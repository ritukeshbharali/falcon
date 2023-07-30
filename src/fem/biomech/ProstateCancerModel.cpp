
/** @file ProstateCancerModel.cpp
 *  @brief Implements a prostate cancer model.
 *  
 *  This class implements a phase-field model for prostate
 *  case (see DOI: 10.1073/pnas.1615791113).
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 24 November 2022
 * 
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022] 
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/geom/Geometries.h>
#include <cstdlib>

/* Include user-defined class headers */

#include "FalconBioMechModels.h"
#include "ProstateCancerModel.h"
#include "util/BasicUtils.h"
#include "util/TbFiller.h"
#include "util/XNames.h"

//=======================================================================
//   class ProstateCancerModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  ProstateCancerModel::SHAPE_PROP                = "shape";
const char*  ProstateCancerModel::DTIME_PROP                = "dtime";
const char*  ProstateCancerModel::LAMBDA_PROP               = "lambda"; 
const char*  ProstateCancerModel::TAU_PROP                  = "tau";
const char*  ProstateCancerModel::GROWTH_RATE_PROP          = "growthRate";
const char*  ProstateCancerModel::APOPTOSIS_RATE_PROP       = "apoptosisRate";
const char*  ProstateCancerModel::EPSILON_PROP              = "epsilon";
const char*  ProstateCancerModel::NUTRIENT_AVG_PROP         = "nutrientAvg";
const char*  ProstateCancerModel::NUTRIENT_DEV_PROP         = "nutrientDev";
const char*  ProstateCancerModel::NUTRIENT_DECAY_PROP       = "nutrientDecay";
const char*  ProstateCancerModel::NUTRIENT_CONSUMPTION_PROP = "nutrientConsumption";
const char*  ProstateCancerModel::GENERATE_NEW_SEED         = "generateNewSeed";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


ProstateCancerModel::ProstateCancerModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );

  const String  context = getContext ();

  // Get the element group assigned to this model.

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();
  rank_   = nodes_ .rank        ();

  // Make sure that the number of spatial dimensions (the rank of the
  // mesh) is valid.

  if ( rank_ < 1 || rank_ > 3 )
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  }

  // Create an internal shape object for computing the element shape
  // functions.

  shape_ = IShapeFactory::newInstance (
    joinNames ( myName_, SHAPE_PROP ),
    conf,
    props
  );

  // Make sure that the rank of the shape matches the rank of the
  // mesh.

  if ( shape_->globalRank() != rank_ )
  {
    throw IllegalInputException (
      context,
      String::format (
        "shape has invalid rank: %d (should be %d)",
        shape_->globalRank (),
        rank_
      )
    );
  }

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements (
    context,
    egroup_.getIndices (),
    shape_->nodeCount  ()
  );


  // Get the DOF space, add displacement and phase-field DOFs

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( 2 );

  // Add dofs

  phiTypes_.ref ( dofTypes_[slice(BEGIN,1)] );
  phiTypes_[0] = dofs_->addType ( "phi" );

  sigTypes_.ref ( dofTypes_[slice(1,END)] );
  sigTypes_[0] = dofs_->addType ( "sig" );
    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  initializeIPMPMap_ ();

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  // Configure parameters

  lambda_       = 1.6e5;
  myProps.find ( lambda_, LAMBDA_PROP );
  myConf. set  ( LAMBDA_PROP, lambda_ );

  tau_          = 0.01;
  myProps.find ( tau_, TAU_PROP );
  myConf. set  ( TAU_PROP, tau_ );

  chi_          = 600.;
  myProps.find ( chi_, GROWTH_RATE_PROP );
  myConf. set  ( GROWTH_RATE_PROP, chi_ );

  A_          = 600.;
  myProps.find ( A_, APOPTOSIS_RATE_PROP );
  myConf. set  ( APOPTOSIS_RATE_PROP, A_ );

  eps_        = 5.0e6;
  myProps.find ( eps_, EPSILON_PROP );
  myConf. set  ( EPSILON_PROP, eps_ );

  savg_       = 1003.75;
  myProps.find ( savg_, NUTRIENT_AVG_PROP );
  myConf. set  ( NUTRIENT_AVG_PROP, savg_ );

  sdev_       = 73.;
  myProps.find ( sdev_, NUTRIENT_DEV_PROP );
  myConf. set  ( NUTRIENT_DEV_PROP, sdev_ );

  delta_       = 1003.75;
  myProps.find ( delta_, NUTRIENT_CONSUMPTION_PROP );
  myConf. set  ( NUTRIENT_CONSUMPTION_PROP, delta_ );

  gamma_       = 1000.;
  myProps.find ( gamma_, NUTRIENT_DECAY_PROP );
  myConf. set  ( NUTRIENT_DECAY_PROP, gamma_ );

  dtime_       = 0.001;
  myProps.find ( dtime_, DTIME_PROP );
  myConf. set  ( DTIME_PROP, dtime_ );

  bool generateNewSeed = false;
  myProps.find ( generateNewSeed, GENERATE_NEW_SEED );

  if ( generateNewSeed )
  {
    // Initialize the seed with time(0), the current unix timestamp
    // On every run, a new set of random numbers are generated.

    srand(time(0));
  }

  // Fill in random values (0.0,1.0]

  srand_.resize( ipCount );
  std::generate(srand_.begin(),srand_.end(), rand );
  srand_ /= max(srand_);

}

ProstateCancerModel::~ProstateCancerModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ProstateCancerModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ProstateCancerModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool ProstateCancerModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  state0;
    Vector  force; 

    // Get the current and old step states.

    StateVector::get       ( state,   dofs_, globdat );
    StateVector::getOld    ( state0,  dofs_, globdat );
    
    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, state, state0 );

    return true;
  }

  /** COMMIT: Requests an action when a computation step has converged.
   *  Ideally, at this point, the material internal state variables
   *  are swapped.
  */ 

  if ( action == Actions::COMMIT )
  {
    return true;
  }

  // CHECK_COMMIT: Can used to discard the current step

  if ( action == Actions::CHECK_COMMIT )
  {
    return true;
  }

  // Post-processing

  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  // Update the stepsize in FE model

  if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void ProstateCancerModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   state,
    const Vector&   state0 )

{

  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // cancer phase-field and nutrient concentration
  
  Vector      phi        ( nodeCount );
  Vector      sig        ( nodeCount );

  // old-step element vector state:
  // cancer phase-field and nutrient concentration
  
  Vector      phi0       ( nodeCount );
  Vector      sig0       ( nodeCount );

  // internal force vector:
  // cancer phase-field and nutrient concentration

  Vector      elemForce1 ( nodeCount );
  Vector      elemForce2 ( nodeCount );

  // element jacobian matrices 
  // (four components)

  Matrix      elemMat1   ( nodeCount, nodeCount ); // phi-phi
  Matrix      elemMat2   ( nodeCount, nodeCount ); // phi-sig
  Matrix      elemMat3   ( nodeCount, nodeCount ); // sig-phi
  Matrix      elemMat4   ( nodeCount, nodeCount ); // sig-sig
  
  // B matrices

  Matrix      be         ( rank_, nodeCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   phiDofs    ( nodeCount );
  IdxVector   sigDofs    ( nodeCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {

    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( phiDofs,  inodes, phiTypes_ );
    dofs_->getDofIndices ( sigDofs,  inodes, sigTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the shape function

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current and old-step element states

    phi  = select ( state,  phiDofs );
    sig  = select ( state,  sigDofs );

    phi0 = select ( state0, phiDofs );
    sig0 = select ( state0, sigDofs );

    // Initialize the internal forces
    // and the jacobian matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0;

    elemMat1   = 0.0;
    elemMat2   = 0.0;
    elemMat3   = 0.0;
    elemMat4   = 0.0;    

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {

      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );

      // Get B-matrix associated with scalar dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the integration point phase-field and nutrient concentration

      Vector Nip           = N   ( ALL, ip   );

      double pf            = dot ( Nip, phi  );
      double sigma         = dot ( Nip, sig  );

      double pf0           = dot ( Nip, phi0 );
      double sigma0        = dot ( Nip, sig0 );

      // Compute derivatives of double well chemical potential

      const double dFdphi    = 32. *  ( pf - 1. ) * pf * ( 2. * pf - 1. );
      const double ddFddphi  = 192. * ( pf * pf - pf ) + 32.;

      // Compute integration point source, based on savg_, smin_, smax_,
      // and the random number srand_[ipoint]

      const double s         = (savg_ - sdev_) + 2.0 * sdev_ * srand_[ipoint]; 

      // Compute jacobian matrix components

      wip         = ipWeights[ip];

      elemMat1   += wip * (  matmul( Nip, Nip ) * ( 1./dtime_ + ddFddphi/tau_ + A_ )
                           + mc2.matmul ( bet, be ) * lambda_ );
      elemMat2   -= wip * (  matmul( Nip, Nip ) * chi_ );
      elemMat3   += wip * (  matmul ( Nip, Nip ) * delta_ ); 
      elemMat4   += wip * (  matmul( Nip, Nip ) * ( 1./dtime_ + gamma_ )
                           + mc2.matmul ( bet, be ) * eps_ );

      // Compute internal forces

      elemForce1 += wip * (  Nip * ( ( pf - pf0 )/dtime_ + dFdphi/tau_ - chi_ * sigma + A_ * pf )
                          +  lambda_ * mc1.matmul( mc2.matmul ( bet, be ), phi ) );
      elemForce2 += wip * (  Nip * ( ( sigma - sigma0 )/dtime_ - s + delta_ * pf + gamma_ * sigma )
                          +  eps_ * mc1.matmul( mc2.matmul ( bet, be ), sig ) );

      
    }  // End of loop on integration points

    // Assembly ...

    mbuilder.addBlock ( phiDofs,  phiDofs, elemMat1 );
    mbuilder.addBlock ( phiDofs,  sigDofs, elemMat2 );
    mbuilder.addBlock ( sigDofs,  phiDofs, elemMat3 );
    mbuilder.addBlock ( sigDofs,  sigDofs, elemMat4 );
    
    select ( force, phiDofs ) += elemForce1;
    select ( force, sigDofs ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool ProstateCancerModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;

  Ref<XTable>  table;
  Vector       weights;
  String       name;


  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table,   ActionParams::TABLE );
  params.get ( name,    ActionParams::TABLE_NAME );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );

  // Element stresses
  if ( name == "source" && 
       table->getRowItems() == elems_.getData() )
  {

    getElemSource_ ( *table, weights);

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getElemSource_
//-----------------------------------------------------------------------


void ProstateCancerModel::getElemSource_

  ( XTable&        table,
    const Vector&  weights )

{

  IdxVector  ielems      = egroup_.getIndices  ();
  const int  ielemCount  = ielems.size         ();
  const int  ipCount     = shape_->ipointCount ();

  Matrix     elSource   ( ielemCount , 1 );

  IdxVector  ielem      ( ielemCount );
  IdxVector  jcols      ( 1   );

  int igpoint = 0;

  // Add the columns for the phase-field to the table.

  jcols[0] = table.addColumn ( "source" );

  elSource      = 0.0;
  ielem         = 0;
  
  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    ielem[ie] = ie;

    for ( int ip = 0; ip < ipCount; ip++, igpoint++ )
    {
      // Store Ip averaged phase-field

      const double s   = (savg_ - sdev_) + 2.0 * sdev_ * srand_[igpoint];
      elSource(ie,0)  += s/ipCount;
    }
  }

  // Add the phase-fields to the table.

  table.addBlock ( ielem, jcols, elSource );

}


//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  ProstateCancerModel::initializeIPMPMap_ ( )

{
  jive::IdxVector   ielems     = egroup_.getIndices  ();

  const idx_t   ielemCount = ielems.size         ();
  const idx_t   ipCount    = shape_->ipointCount ();

        idx_t   ipoint     = 0;

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // loop over integration points 

    for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
    {
      ipMpMap_ ( ielem, ip ) = ipoint;
    }
  }
}


//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void ProstateCancerModel::setStepSize_

  ( const Properties&  params )

{
  double       dt;
  params.get ( dt,  XProps::STEP_SIZE   );

  System::debug() << "Setting step size to " << dt << " ("
    << dt / dtime_ * 100 << "\% of previous step size)" << "\n";

  dtime_ = dt;
  
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newProstateCancerModel
//-----------------------------------------------------------------------


static Ref<Model>     newProstateCancerModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<ProstateCancerModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareProstateCancerModel
//-----------------------------------------------------------------------


void declareProstateCancerModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "ProstateCancer", & newProstateCancerModel );
}
