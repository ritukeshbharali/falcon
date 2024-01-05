
/** @file CahnHilliardModel.cpp
 *  @brief Implements a Cahn_Hilliard model.
 *  
 *  This class implements a Cahn-Hilliard model
 *  using operator split.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 23 October 2023
 * 
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022] 
 * 
 *  To-do: Replace double well chemical potential
 *         with a user-defined option.
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/geom/Geometries.h>
#include <cstdlib>

/* Include user-defined class headers */

#include "FalconBasicModels.h"
#include "CahnHilliardModel.h"
#include "util/XNames.h"

//=======================================================================
//   class CahnHilliardModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  CahnHilliardModel::SHAPE_PROP                = "shape";
const char*  CahnHilliardModel::DTIME_PROP                = "dtime";
const char*  CahnHilliardModel::DIFFUSIVITY_PROP          = "diffusivity";
const char*  CahnHilliardModel::GAMMA_PROP                = "gamma";
const char*  CahnHilliardModel::GENERATE_NEW_SEED         = "generateNewSeed";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


CahnHilliardModel::CahnHilliardModel

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


  // Get the DOF space

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( 2 );

  // Add concentration and chemical potential dofs

  concTypes_.ref ( dofTypes_[slice(BEGIN,1)] );
  concTypes_[0] = dofs_->addType ( "c" );

  cPotTypes_.ref ( dofTypes_[slice(1,END)] );
  cPotTypes_[0] = dofs_->addType ( "mu" );
    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  // Configure parameters

  dtime_       = 5.e-6;
  myProps.find ( dtime_, DTIME_PROP );
  myConf. set  ( DTIME_PROP, dtime_ );

  D_            = 1.0;
  myProps.find ( D_, DIFFUSIVITY_PROP );
  myConf. set  ( DIFFUSIVITY_PROP, D_ );

  gamma_       = 0.01;
  myProps.find ( gamma_, GAMMA_PROP );
  myConf. set  ( GAMMA_PROP, gamma_ );

  bool generateNewSeed = false;
  myProps.find ( generateNewSeed, GENERATE_NEW_SEED );

  if ( generateNewSeed )
  {
    // Initialize the seed with time(0), the current unix timestamp
    // On every run, a new set of random numbers are generated.

    srand(time(0));
  }

}

CahnHilliardModel::~CahnHilliardModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void CahnHilliardModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void CahnHilliardModel::getConfig ( const Properties& conf,
                                    const Properties& globdat )  const
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool CahnHilliardModel::takeAction

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


void CahnHilliardModel::getMatrix_

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
  // concentration and chemical potential
  
  Vector      conc        ( nodeCount );
  Vector      cPot        ( nodeCount );

  // old-step element vector state:
  // concentration
  
  Vector      conc0       ( nodeCount );

  // internal force vector:
  // concentration and chemical potential

  Vector      elemForce1 ( nodeCount );
  Vector      elemForce2 ( nodeCount );

  // element jacobian matrices 
  // (four components)

  Matrix      elemMat1   ( nodeCount, nodeCount ); // conc-conc
  Matrix      elemMat2   ( nodeCount, nodeCount ); // conc-cPot
  Matrix      elemMat3   ( nodeCount, nodeCount ); // cPot-conc
  Matrix      elemMat4   ( nodeCount, nodeCount ); // cPot-cPot
  
  // B matrices

  Matrix      be         ( rank_, nodeCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   concDofs   ( nodeCount );
  IdxVector   cPotDofs   ( nodeCount );

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
    dofs_->getDofIndices ( concDofs, inodes, concTypes_ );
    dofs_->getDofIndices ( cPotDofs, inodes, cPotTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the shape function

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current and old-step element states

    conc  = select ( state,  concDofs );
    cPot  = select ( state,  cPotDofs );

    conc0 = select ( state0, concDofs );

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

      // Get B-matrix associated with scalar dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the integration point concentration and chemical potential

      Vector Nip  = N   ( ALL, ip   );

      const double c    = dot ( Nip, conc  );
      const double mu   = dot ( Nip, cPot  );
      const double c0   = dot ( Nip, conc0 );

      // Compute derivatives of double well chemical potential

      // const double F    = 100. * ::pow( c,2.0 ) * ::pow( 1.0-c,2.0 );
      const double dFdc    = 200. * ( c - 1. ) * c * ( 2. * c - 1. );
      const double ddFddc  = 1200.* ( c * c - c ) + 200.; 

      // Compute jacobian matrix components

      wip         = ipWeights[ip];

      elemMat1   += wip / dtime_ * matmul( Nip, Nip );
      elemMat2   += wip * D_ * mc2.matmul ( bet, be );
      elemMat3   -= wip * (  matmul ( Nip, Nip ) * ddFddc
                           - gamma_ * mc2.matmul ( bet, be ) ); 
      elemMat4   += wip * matmul( Nip, Nip );

      // Compute internal forces

      elemForce1 += wip * (  Nip * ( ( c - c0 )/dtime_ )
                          +  D_ * mc1.matmul( mc2.matmul ( bet, be ), cPot ) );
      elemForce2 += wip * (  Nip * ( mu - dFdc )
                          +  gamma_ * mc1.matmul( mc2.matmul ( bet, be ), conc ) );

      
    }  // End of loop on integration points

    // Assembly ...

    mbuilder.addBlock ( concDofs,  concDofs, elemMat1 );
    mbuilder.addBlock ( concDofs,  cPotDofs, elemMat2 );
    mbuilder.addBlock ( cPotDofs,  concDofs, elemMat3 );
    mbuilder.addBlock ( cPotDofs,  cPotDofs, elemMat4 );
    
    select ( force, concDofs ) += elemForce1;
    select ( force, cPotDofs ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void CahnHilliardModel::setStepSize_

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
//   newCahnHilliardModel
//-----------------------------------------------------------------------


static Ref<Model>     newCahnHilliardModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<CahnHilliardModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareCahnHilliardModel
//-----------------------------------------------------------------------


void declareCahnHilliardModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "CahnHilliard", & newCahnHilliardModel );
}
