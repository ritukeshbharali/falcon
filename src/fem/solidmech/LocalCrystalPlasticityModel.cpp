
/** @file LocalCrystalPlasticityModel.cpp
 *  @brief Crystal visco-plasticity model without regularization.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 04 June 2024
 *
 *  Updates (when, what and who)
 *     - [14 June 2024] added functions to write stress,
 *       strain, and slip nodal and element tables. (RB)
 *     - [25 June 2024] added user-defined rotation 
 *       parameter operating on the slip plane and 
 *       direction for flexibility. (RB)
 * 
 *  @todo Simplify some data structures, reduce loops.
 *
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/Float.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/implict/ArclenActions.h>

/* Include other headers */

#include "FalconSolidMechModels.h"
#include "LocalCrystalPlasticityModel.h"
#include "util/TbFiller.h"
#include "util/XNames.h"
#include "util/Constants.h"
#include "util/MathUtils.h"
#include "util/TensorUtils.h"

#include "materials/Material.h"
#include "materials/HookeMaterial.h"


//=======================================================================
//   class LocalCrystalPlasticityModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  LocalCrystalPlasticityModel::DISP_NAMES[3]  = { "dx", "dy", "dz" };

const char*  LocalCrystalPlasticityModel::SHAPE_PROP     = "shape";
const char*  LocalCrystalPlasticityModel::MATERIAL_PROP  = "material";
const char*  LocalCrystalPlasticityModel::DTIME_PROP     = "dtime";
const char*  LocalCrystalPlasticityModel::SLIPS_PROP     = "slips";
const char*  LocalCrystalPlasticityModel::ROTATION_PROP  = "rotation";
const char*  LocalCrystalPlasticityModel::RHO_PROP       = "rho";
const char*  LocalCrystalPlasticityModel::IPNODES_PROP   = "ipNodes";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


LocalCrystalPlasticityModel::LocalCrystalPlasticityModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super   ( name ),
    nslips_ ( 2 ),
    rot_    ( 0.0 ),
    dtime_  ( 1.0 ),
    tstar_  ( 1.0 ),
    tauY_   ( 1.0 ),
    n_      ( 1.0 )

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

  // Get the number of slip systems for the model

  myProps.get ( slips_, SLIPS_PROP );
  myConf. set ( SLIPS_PROP, slips_ );

  nslips_ = slips_.size();

  if ( nslips_ < 1 )
  {
    throw IllegalInputException (
      context,
      "at least one slip system must exist!"
    );
  }

  // Get the DofSpace, add displacement, slip, tau dof types
  // Displacement and slip operates on nodes
  // tau operates on ip nodes (special!)

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ + nslips_ );

  // Make a shallow copy of the first (displacement) part of 
  // dofTypes_ and add displacement dof names

  dispTypes_.ref ( dofTypes_[slice(BEGIN,rank_)] );

  for ( int i = 0; i < rank_; i++ )
  {
    dispTypes_[i] = dofs_->addType ( DISP_NAMES[i] );
  }

  // Make a shallow copy of dofTypes_ for slip dofs

  slipTypes_.ref ( dofTypes_[slice(rank_, rank_ + nslips_)] );

  for ( int i = 0; i < nslips_; i++ )
  {
    String slipName = "slip" + String(i);
    slipTypes_[i] = dofs_->addType ( slipName );
  }

  // Assign displacement dofs to each node
  
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_[slice(0, rank_)]
  );

  // Get integration point (dummy) nodes

  myProps.get ( ipNGroup_, IPNODES_PROP );

  ipnodes_  = NodeGroup::get ( ipNGroup_, nodes_, globdat, context );

  // Compute the total number of integration points.

  const int ipCount = shape_->ipointCount() * egroup_.size();

  if ( ipnodes_.size() != ipCount )
  {
    throw IllegalInputException (
      context,
      "wrong ipNodes size!"
    );
  }

  myConf.set ( IPNODES_PROP, ipNGroup_ );

  dofs_->addDofs (
    ipnodes_.getIndices(),
    dofTypes_[slice(rank_, END)]
  );

  // Do we need to store strains? Or compute them on-the-fly?
  // Store all tau

  strain_.resize ( STRAIN_COUNTS[rank_], ipCount );
  strain_ = 0.0;

  // Create a material model object. Ensure it is a Hooke material 
  // with dynamicCast

  Ref<Material> tmpMat = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  material_            = dynamicCast<HookeMaterial> ( tmpMat );

  if ( material_ == NIL )
  {
    throw IllegalInputException (
      context,
      "material type is not a Hooke material"
    );
  }

  initializeIPMPMap_ ();

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  // Set up the variables for the slip systems

  sm_   .resize ( nslips_ );
  Dsm_  .resize ( nslips_ );
  Esm_  .resize ( nslips_, nslips_ );

  tstar_.resize ( nslips_ );
  tauY_ .resize ( nslips_ );
  n_    .resize ( nslips_ );

  // Allocate for slip - direction dyadics
  
  for ( int islip = 0; islip < nslips_; islip++ )
  {
    sm_[islip].resize( STRAIN_COUNTS[rank_] );
  }

}


LocalCrystalPlasticityModel::~LocalCrystalPlasticityModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  MChain1     mc1;

  const String  context = getContext ();

  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  myProps.get ( dtime_, DTIME_PROP );

  material_->configure ( matProps, globdat );

  D0_.resize( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  D0_ = material_-> getStiffMat();

  // Find rotation angle (in degrees) and convert to radian
  
  myProps.find ( rot_, ROTATION_PROP ); rot_ *= PI/180;

  // Create a rotation matrix

  Matrix R (3,3); R = 0.0;
  R(0,0) = R(1,1) = cos(rot_);
  R(0,1) = -sin(rot_);
  R(1,0) =  sin(rot_);
  R(2,2) =  1.0;

  for ( int i = 0; i < nslips_; i++ )
  {
    String slipName = slips_[i];

    Properties mySlipProps = myProps.getProps ( slipName );

    mySlipProps.get ( tstar_[i], "tstar" );
    mySlipProps.get ( tauY_[i],  "tauY" );
    mySlipProps.get ( n_[i],     "n" );

    Vector plane ( 3 );    // 3D
    Vector dir   ( 3 );    // 3D

    mySlipProps.get ( plane, "plane" );
    mySlipProps.get ( dir,   "direction" );

    if ( plane.size() != 3 )
    {
      throw IllegalInputException (
      context,
      "plane must be a vector size 3"
      );
    }

    if ( dir.size() != 3 )
    {
      throw IllegalInputException (
      context,
      "direction must be a vector size 3"
      );
    }

    // Rotate slip plane and direction

    plane = mc1.matmul( R, plane );
    dir   = mc1.matmul( R, dir   );

    // Normalize plane and direction

    double fac1 = 1.0/jem::numeric::norm2( plane );
    double fac2 = 1.0/jem::numeric::norm2( dir );

    plane *= fac1;
    dir   *= fac2;

    // Compute sm_

    Matrix tmp = tensorUtils::otimes( plane, dir );
    sm_[i]     = tensorUtils::tensor2voigtStrain( tmp, 
                                STRAIN_COUNTS[rank_] );

  }

  // Compute some more coefficients

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    Dsm_[islip].resize( STRAIN_COUNTS[rank_] );

    matmul ( Dsm_[islip], D0_, sm_[islip] );

    for ( int jslip = 0; jslip < nslips_; jslip++ )
    {
      Esm_(islip,jslip) = dot( Dsm_[islip], sm_[jslip] );
    }
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP );

  myConf. set ( DTIME_PROP, dtime_ );

  for ( int i = 0; i < nslips_; i++ )
  {
    String slipName = slips_[i];

    Properties mySlipConf = myConf.makeProps ( slipName );

    mySlipConf.set ( "tstar", tstar_[i] );
    mySlipConf.set ( "tauY",  tauY_[i] );
    mySlipConf.set ( "n",     n_[i] );
  }

  material_->getConfig ( matConf, globdat );

  myConf. set ( ROTATION_PROP, rot_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool LocalCrystalPlasticityModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;
  using jive::implict::ArclenActions;

  // Compute the internal force vector

  if ( action == Actions::GET_INT_VECTOR )
  {

    Vector  state, state0;
    Vector  force; 

    // Get the current and previous step state.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );
    
    // Get the internal force vector.

    params.get ( force, ActionParams::INT_VECTOR );

    getMatrix_ ( nullptr, force, state, state0 );

    return true;
  }

  // Compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  state, state0;
    Vector  force; 

    // Get the current and previous step state.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );
    
    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( mbuilder, force, state, state0 );

    return true;
  }

  // Compute mass matrix 

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );
    
    getMatrix2_( *mbuilder );

    return true;
  }

  // Perform a commit once the solver converges and the solution
  // is accepted. Ideally, at this point, state variables are swapped.

  if ( action == Actions::COMMIT )
  {
    return true;
  }

  // CHECK_COMMIT: Can used to discard the current step 

  if ( action == Actions::CHECK_COMMIT )
  {
    //checkCommit_ ( params );
    return true;
  }

  /** GET_TABLE: Requests a post-processing action, such as computing
   *  tabular data (stress, strains, etc.) to print VTK files.
  */ 
  
  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       state,
    const Vector&       state0 )

{
  // Get the elements associated with this model
  IdxVector   ielems     = egroup_ .getIndices ();

  // Get the integration point nodes for this model
  IdxVector   ipnodes    = ipnodes_.getIndices ();

  // Get some additional details
  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];
  const int   dispCount  = nodeCount * rank_;

  // Allocate grads and coords
  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // Allocate vectors to store element states
  // displacement (conventional nodal qtys)
  Vector         disp   ( dispCount );

  // Allocate vector for slip dofs. These are a bit 
  // tricky. For each slip system we extract a
  // ip nodal vector (i.e, inner vector has size
  // same as ipCount)
  vector<Vector> slip   ( nslips_ );
  vector<Vector> slip0  ( nslips_ );

  Vector         tau    ( nslips_ );

  // Allocate matrix and vectors for material update
  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );

  // Allocate vectors to store element internal
  // forces
  Vector         elemForce1 ( dispCount );  // disp
  vector<Vector> elemForce2 ( nslips_   );  // slip

  // element stiffness matrices 
  // (4 block components )

  // Displacement row

  Matrix          elemMat11  ( dispCount, dispCount );
  vector<Matrix>  elemMat12  ( nslips_ );

  // Slip rows

  vector<Matrix>  elemMat21  ( nslips_ );
  vector<Matrix>  elemMat22  ( nslips_ );

  // resize inner vectors and matrices
  // i.e., corresponding to each slip

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    // Element slips

    slip[islip] .resize( ipCount );
    slip0[islip].resize( ipCount );

    // Element slip internal force

    elemForce2[islip].resize( ipCount );

    // Element slip and tau matrices

    elemMat12[islip].resize ( dispCount, ipCount );

    elemMat21[islip].resize ( ipCount, dispCount );
    elemMat22[islip].resize ( ipCount, ipCount   );
  }
  
  // B matrices
  // displacement
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  // Element nodes and corresponding dofs
  
  IdxVector         inodes   ( nodeCount );
  IdxVector         dispDofs ( dispCount );
  vector<IdxVector> slipDofs ( nslips_ );

  // re-size the inner vectors

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    slipDofs[islip].resize( ipCount );
  }

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

    const int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );
    
    // Get current element displacements

    disp   = select ( state,  dispDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix to zero
    
    elemForce1 = 0.0;

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      elemForce2[islip] = 0.0;
    }

    if ( mbuilder != nullptr )
    {
      elemMat11   = 0.0;

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        elemMat12[islip]   = 0.0;

        elemMat21[islip]   = 0.0;
        elemMat22[islip]   = 0.0;
      }    
    }

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      const int ipoint = ipMpMap_ ( ielem, ip );

      // We use this information later, ipnodes[ipoint] is the
      // dummy integration point node, where tau is defined.

      // Compute the displacement B-matrix for this integration point.

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // Compute the strain for this integration point.

      matmul ( strain,  bd, disp  );

      // Store the regular strain components.

      strain_(slice(BEGIN,strCount),ipoint) = strain;

      // Compute the elastic stress

      material_->update ( stress, stiff, strain, 0  );

      // Compute the slips for this ipoint

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        slipDofs[islip][ip] = dofs_-> getDofIndex ( 
                  ipnodes[ipoint], slipTypes_[islip] );

        slip[islip]   = select ( state,  slipDofs[islip] );
        slip0[islip]  = select ( state0, slipDofs[islip] );

        // Update stress

        stress -= slip[islip][ip] * Dsm_[islip];

        // Compute tau

        tau[islip] = dot( sm_[islip], stress );
      }

      // Compute the overstress function Phi and its derivative
      // Assume n = 1 for a start!

      Vector Phi  ( nslips_ );
      Vector dPhi ( nslips_ );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        Phi[islip]  = sign( tau[islip] ) / tstar_[islip]
                     * ::pow( fabs( tau[islip] ) / 
                              tauY_[islip], n_[islip] );

        dPhi[islip] = 1.0 / tstar_[islip]
                        * n_[islip] / tauY_[islip]
                        * ::pow( fabs( tau[islip] ) / 
                                  tauY_[islip], n_[islip] - 1.0 );

        if ( fabs(dPhi[islip]) < 1.e-15 )
        {
          // Provide a dummy value
          dPhi[islip] = 1.e-10;
        }
      }

      // Compute weight of ip

      wip         = ipWeights[ip];

      // -----------------------------------------
      // Compute the internal force components
      // -----------------------------------------

      elemForce1 += wip * mc1.matmul( bdt, stress );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        elemForce2[islip]  = wip * ( ( slip[islip][ip] - slip0[islip][ip] ) / dtime_
                                   - Phi[islip]  
                                   );

        select ( force, slipDofs[islip] ) = elemForce2[islip];

        // Note: Esm_ cross hardening is not considered here!
      }

      // -----------------------------------------
      // Compute the stiffness matrix components
      // -----------------------------------------

      // Displacement row

      elemMat11  += wip * mc3.matmul ( bdt, stiff, bd );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        // Displacement row

        elemMat12[islip](ALL,ip) -= wip * mc1.matmul ( bdt, Dsm_[islip] );

        // Slip row

        elemMat21[islip](ip,ALL)  = -wip * dPhi[islip]/dtime_ * mc1.matmul ( Dsm_[islip], bd );
        elemMat22[islip](ip,ip)   =  wip * ( 1.0 + dPhi[islip] * Esm_(islip,islip) );

        if ( mbuilder != nullptr )
        {
          mbuilder -> setBlock ( slipDofs[islip], dispDofs,         elemMat21[islip] );
          mbuilder -> setBlock ( slipDofs[islip], slipDofs[islip],  elemMat22[islip] );
        }
      }

    }  // End of loop on integration points

    // Assembly ...

    if ( mbuilder != nullptr )
    {
      mbuilder -> addBlock ( dispDofs, dispDofs, elemMat11 );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        mbuilder -> addBlock ( dispDofs, slipDofs[islip], elemMat12[islip] );
      }      
    }

    select ( force, dispDofs ) += elemForce1;
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void LocalCrystalPlasticityModel::getMatrix2_

    ( MatrixBuilder&          mbuilder )
{
  IdxVector   ielems     = egroup_.getIndices  ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   dofCount   = rank_ * nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );
  
  Matrix      elemMat    ( dofCount, dofCount );

  Matrix      R          ( rank_, rank_ );

  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount );
  Matrix      Nt         = N.transpose ( ); 

  IdxVector   inodes     ( nodeCount );
  IdxVector   idofs      ( dofCount  );

  Vector      ipWeights  ( ipCount   );

  MChain3     mc3;

  double      rho = 0.0;

  R = 0.0;
 
  for ( int i = 0; i < rank_ ; i++ )
  {
    R(i,i) = rho;
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dispTypes_ );    // dofTypes_

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Assemble the element matrix and the internal force vector.

    elemMat   = 0.0;

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // compute matrix of shpae function N       

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Add the element secant matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}

//-----------------------------------------------------------------------
//   getArcFunc_
//-----------------------------------------------------------------------

void LocalCrystalPlasticityModel::getArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  LocalCrystalPlasticityModel::initializeIPMPMap_ ( )

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
//   getTable_
//-----------------------------------------------------------------------


bool LocalCrystalPlasticityModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  Ref<XTable>  table;
  Vector       weights;
  String       name;


  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table,   ActionParams::TABLE );
  params.get ( name,    ActionParams::TABLE_NAME );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );

  // Nodal stresses
  if ( name == "stress" && 
       table->getRowItems() == nodes_.getData() )
  {
    Vector  state;
    StateVector::get ( state,     dofs_, globdat  );

    getStress_ ( *table, weights, state);

    return true;
  }

  // Element stresses
  if ( name == "stress" && 
       table->getRowItems() == elems_.getData() )
  {
    Vector  state;
    StateVector::get ( state,     dofs_, globdat  );

    getElemStress_ ( *table, weights, state);

    return true;
  }

  // Nodal strains
  if ( name == "strain" && 
       table->getRowItems() == nodes_.getData() )
  {
    getStrain_ ( *table, weights);

    return true;
  }

  // Element strains
  if ( name == "strain" && 
       table->getRowItems() == nodes_.getData() )
  {
    getElemStrain_ ( *table, weights);

    return true;
  }

  // Nodal strains
  if ( name == "slip" && 
       table->getRowItems() == nodes_.getData() )
  {
    Vector  state;
    StateVector::get ( state, dofs_, globdat  );

    getSlip_ ( *table, weights, state);

    return true;
  }

  // Nodal strains
  if ( name == "slip" && 
       table->getRowItems() == elems_.getData() )
  {
    Vector  state;
    StateVector::get ( state, dofs_, globdat  );

    getElemSlip_ ( *table, weights, state);

    return true;
  }

  // 
  if ( name == "xoutTable" )
  {
    Vector  state;
    String  contents;

    StateVector::get ( state,     dofs_, globdat  );
    params.      get ( contents, "contentString" );

    
    getOutputData_ ( table, weights, contents, state );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getStress_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  state )

{
  IdxVector  ielems     = egroup_ .getIndices  ();
  IdxVector  ipnodes    = ipnodes_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];

  // Get the shape function

  Matrix N  = shape_->getShapeFunctions ();

  Matrix     ndStress   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  Matrix     stiff      ( strCount, strCount );
  Vector     stressIp   ( strCount );
  vector<Vector> slip   ( nslips_ );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  vector<IdxVector> slipDofs ( nslips_ );

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    slip[islip]     .resize( ipCount );
    slipDofs[islip] .resize( ipCount );
  }

  MChain1    mc1;

  // Add the columns for the stress components to the table.

  switch ( strCount )
  {
  case 1:

    jcols[0] = table.addColumn ( "stress_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_xy" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );

    break;  

  case 6:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );
    jcols[4] = table.addColumn ( "stress_yz" );
    jcols[5] = table.addColumn ( "stress_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount )
    );
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    const int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      const int ipoint = ipMpMap_ ( ielem, ip );

      // Extrapolate the integration point stresses to the nodes using
      // the transposed shape functions.

      material_->update ( stressIp, stiff, strain_(ALL,ipoint), 0  );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        slipDofs[islip][ip] = dofs_-> getDofIndex ( 
                  ipnodes[ipoint], slipTypes_[islip] );

        slip[islip]   = select ( state,  slipDofs[islip] );

        // Update stress

        stressIp -= slip[islip][ip] * Dsm_[islip];
      }

      ndStress  += matmul ( N(ALL,ip), stressIp );
      ndWeights += N(ALL,ip);
    }

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols, ndStress );
  }
}


//-----------------------------------------------------------------------
//   getElemStress_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getElemStress_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  state )

{
  IdxVector  ielems     = egroup_ .getIndices  ();
  IdxVector  ipnodes    = ipnodes_.getIndices  ();

  const int  elemCount   = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];

  // Get the shape function

  Matrix N  = shape_->getShapeFunctions ();

  Matrix     elStress   ( elemCount, strCount );

  Matrix     stiff      ( strCount, strCount );
  Vector     stressIp   ( strCount );
  vector<Vector> slip   ( nslips_ );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  vector<IdxVector> slipDofs ( nslips_ );

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    slip[islip]     .resize( ipCount );
    slipDofs[islip] .resize( ipCount );
  }

  MChain1    mc1;

  // Add the columns for the stress components to the table.

  switch ( strCount )
  {
  case 1:

    jcols[0] = table.addColumn ( "stress_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_xy" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );

    break;  

  case 6:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );
    jcols[4] = table.addColumn ( "stress_yz" );
    jcols[5] = table.addColumn ( "stress_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount )
    );
  }

  elStress = 0.0;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    const int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      const int ipoint = ipMpMap_ ( ielem, ip );

      // Extrapolate the integration point stresses to the nodes using
      // the transposed shape functions.

      material_->update ( stressIp, stiff, strain_(ALL,ipoint), 0  );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        slipDofs[islip][ip] = dofs_-> getDofIndex ( 
                  ipnodes[ipoint], slipTypes_[islip] );

        slip[islip]   = select ( state,  slipDofs[islip] );

        // Update stress

        stressIp -= slip[islip][ip] * Dsm_[islip];
      }

      elStress(ie,ALL) += stressIp/ipCount;
    }
  }

  // Add the stresses to the table.

  table.addBlock ( ielems, jcols, elStress );
}


//-----------------------------------------------------------------------
//   getStrain_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getStrain_

  ( XTable&        table,
    const Vector&  weights )

{
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount = ielems.size         ();
  const int  nodeCount  = shape_->nodeCount   ();
  const int  ipCount    = shape_->ipointCount ();
  const int  strCount   = STRAIN_COUNTS[rank_];

  Matrix     N          = shape_->getShapeFunctions ();

  Matrix     ndStrain   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  // Add the columns for the normal strain components to the table.

  switch ( strCount )
  {
  case 1:

    jcols[0] = table.addColumn ( "e_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_xy" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_zz" );
    jcols[3] = table.addColumn ( "e_xy" );

    break;  

  case 6:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_zz" );
    jcols[3] = table.addColumn ( "e_xy" );
    jcols[4] = table.addColumn ( "e_yz" );
    jcols[5] = table.addColumn ( "e_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of strain components: " +
      String ( strCount )
    );
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    // Get the element nodes.

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndStrain  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      const int ipoint = ipMpMap_ ( ielem, ip );

      // Extrapolate the integration point strains to the nodes using
      // the transposed shape functions.

      ndStrain  += matmul ( N(ALL,ip), strain_(slice(BEGIN,strCount),ipoint) );
      ndWeights += N(ALL,ip);
    }

    // Increment the table weights. When the complete table has been
    // filled, Jive will divide each row in the table by the
    // corresponding table weight. In this way the strain components
    // are automatically averaged over all elements that are attached
    // to a node. The weight vector is initially zero.

    select ( weights, inodes ) += ndWeights;

    // Add the strains to the table.

    table.addBlock ( inodes, jcols, ndStrain );
  }
}


//-----------------------------------------------------------------------
//   getElemStrain_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getElemStrain_

  ( XTable&        table,
    const Vector&  weights )

{
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  elemCount  = ielems.size         ();
  const int  nodeCount  = shape_->nodeCount   ();
  const int  ipCount    = shape_->ipointCount ();
  const int  strCount   = STRAIN_COUNTS[rank_];

  Matrix     elStrain   ( elemCount, strCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  // Add the columns for the normal strain components to the table.

  switch ( strCount )
  {
  case 1:

    jcols[0] = table.addColumn ( "e_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_xy" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_zz" );
    jcols[3] = table.addColumn ( "e_xy" );

    break;  

  case 6:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_zz" );
    jcols[3] = table.addColumn ( "e_xy" );
    jcols[4] = table.addColumn ( "e_yz" );
    jcols[5] = table.addColumn ( "e_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of strain components: " +
      String ( strCount )
    );
  }

  elStrain = 0.0;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    // Get the element nodes.

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      const int ipoint = ipMpMap_ ( ielem, ip );

      // Compute integration point-averaged strains

      elStrain(ie,ALL) += strain_(slice(BEGIN,strCount),ipoint)/ipCount;
    }
  }

  // Add the strains to the table.

  table.addBlock ( ielems, jcols, elStrain );
}


//-----------------------------------------------------------------------
//   getSlip_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getSlip_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  state )

{
  IdxVector  ielems     = egroup_ .getIndices  ();
  IdxVector  ipnodes    = ipnodes_.getIndices  ();

  const int  elemCount  = ielems.size         ();
  const int  nodeCount  = shape_->nodeCount   ();
  const int  ipCount    = shape_->ipointCount ();

  Matrix     N          = shape_->getShapeFunctions ();

  Matrix     ndSlip     ( nodeCount, nslips_ );
  Vector     ndWeights  ( nodeCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( nslips_   );

  // Add the columns to the table.

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    String tauName = "slip" + String(islip);
    jcols[islip]   = table.addColumn ( tauName );
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    // Get the element nodes.

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndSlip    = 0.0;
    ndWeights = 0.0;

    IdxVector slipDofs ( nslips_ );
    Vector    slip     ( nslips_ );

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );
      
      dofs_     -> getDofIndices ( slipDofs, ipnodes[ipoint], slipTypes_ );
      slip       = state[slipDofs];

      // Extrapolate the integration point taus to the nodes using
      // the transposed shape functions.

      ndSlip    += matmul ( N(ALL,ip), slip );
      ndWeights += N(ALL,ip);
    }

    // Increment the table weights. When the complete table has been
    // filled, Jive will divide each row in the table by the
    // corresponding table weight. In this way the strain components
    // are automatically averaged over all elements that are attached
    // to a node. The weight vector is initially zero.

    select ( weights, inodes ) += ndWeights;

    // Add the strains to the table.

    table.addBlock ( inodes, jcols, ndSlip);
  }
}


//-----------------------------------------------------------------------
//   getElemSlip_
//-----------------------------------------------------------------------


void LocalCrystalPlasticityModel::getElemSlip_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  state )

{
  IdxVector  ielems     = egroup_ .getIndices  ();
  IdxVector  ipnodes    = ipnodes_.getIndices  ();

  const int  elemCount  = ielems.size         ();
  const int  nodeCount  = shape_->nodeCount   ();
  const int  ipCount    = shape_->ipointCount ();

  Matrix     elSlip     ( elemCount, nslips_ );
  Vector     ndWeights  ( nodeCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( nslips_   );

  // Add the columns to the table.

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    String slipName = "slip" + String(islip);
    jcols[islip]   = table.addColumn ( slipName );
  }

  elSlip = 0.0;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    // Get the element nodes.

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    IdxVector slipDofs ( nslips_ );
    Vector    slip     ( nslips_ );

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      const int ipoint = ipMpMap_ ( ielem, ip );
      
      dofs_     -> getDofIndices ( slipDofs, ipnodes[ipoint], slipTypes_ );
      slip       = state[slipDofs];

      // Compute integration point-averaged slip

      elSlip(ie,ALL) +=  slip/ipCount;
    }
  }

  // Add the strains to the table.

  table.addBlock ( ielems, jcols, elSlip);
}


//-----------------------------------------------------------------------
//   getOutputData_
//-----------------------------------------------------------------------

/** all data at Gauss points are computed, however, only those requested
 *  in the vtk block in the .pro file will be written to vtu files.
 *  Example: vtk.data = "stress_yy | stress_xx" for xx and yy stress 
 *  components.
 */

void LocalCrystalPlasticityModel::getOutputData_

  ( Ref<XTable>    table,
    const Vector&  weights,
    const String&  contents,
    const Vector&  state )
          
{
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void LocalCrystalPlasticityModel::checkCommit_

  ( const Properties&  params )

{
  // System::info() << myName_ << " : check commit ... do nothing!\n";
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newLocalCrystalPlasticityModel
//-----------------------------------------------------------------------


static Ref<Model>     newLocalCrystalPlasticityModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<LocalCrystalPlasticityModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareLocalCrystalPlasticityModel
//-----------------------------------------------------------------------


void declareLocalCrystalPlasticityModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "LocalCrystalPlasticity", & newLocalCrystalPlasticityModel );
}
