
/** @file CrystalViscoPlasticityModel.cpp
 *  @brief Crystal visco-plasticity model with gradient regularization.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 08 April 2024
 *
 *  Updates (when, what and who)
 * 
 *  TO-DO: Performance enhancements! Too many loops!
 *  Initialize variables together for better cache
 *  utilization!
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
#include "CrystalViscoPlasticityModel.h"
#include "util/TbFiller.h"
#include "util/XNames.h"
#include "util/Constants.h"
#include "util/MathUtils.h"
#include "util/TensorUtils.h"

#include "materials/Material.h"
#include "materials/HookeMaterial.h"


//=======================================================================
//   class CrystalViscoPlasticityModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  CrystalViscoPlasticityModel::DISP_NAMES[3]  = { "dx", "dy", "dz" };

const char*  CrystalViscoPlasticityModel::SHAPE_PROP     = "shape";
const char*  CrystalViscoPlasticityModel::MATERIAL_PROP  = "material";
const char*  CrystalViscoPlasticityModel::DTIME_PROP     = "dtime";
const char*  CrystalViscoPlasticityModel::SLIPS_PROP     = "slips";
const char*  CrystalViscoPlasticityModel::RHO_PROP       = "rho";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


CrystalViscoPlasticityModel::CrystalViscoPlasticityModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super   ( name ),
    nslips_ ( 2 ),
    dtime_  ( 1.0 ),
    tstar_  ( 1.0 ),
    tauY_   ( 1.0 ),
    n_      ( 1.0 )

{
  System::out() << "Enters constructor \n";

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

  // Get the DofSpace, add displacement and slip dof types

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ + 2 * nslips_ );

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

  // Make a shallow copy of dofTypes_ for tau dofs

  tauTypes_.ref ( dofTypes_[slice(rank_ + nslips_, END)] );

  for ( int i = 0; i < nslips_; i++ )
  {
    String tauName = "tau" + String(i);
    tauTypes_[i] = dofs_->addType ( tauName );
  }

  // Assign dofs to each node
  
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  // Do we need to store strains? Or compute them on-the-fly?

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
  
  for ( int islip = 0; islip < nslips_; islip++ )
  {
    sm_[islip].resize( STRAIN_COUNTS[rank_] );
  }  
}


CrystalViscoPlasticityModel::~CrystalViscoPlasticityModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void CrystalViscoPlasticityModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  const String  context = getContext ();

  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  myProps.get ( dtime_, DTIME_PROP );

  material_->configure ( matProps, globdat );

  D0_.resize( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  D0_ = material_-> getStiffMat();

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


void CrystalViscoPlasticityModel::getConfig ( const Properties& conf,
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
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool CrystalViscoPlasticityModel::takeAction

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

  // compute the tangent stiffness matrix and internal force vector

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

  // compute mass matrix 

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );
    
    getMatrix2_( *mbuilder );

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

  /** CHECK_COMMIT: Can used to discard the current step
  */ 

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


void CrystalViscoPlasticityModel::getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       state,
    const Vector&       state0 )

{
  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of displacement dofs

  const int   dispCount  = nodeCount * rank_;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement, slips and tau
  
  Vector         disp   ( dispCount );
  vector<Vector> slip   ( nslips_ );
  vector<Vector> slip0  ( nslips_ );
  vector<Vector> tau    ( nslips_ );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      stressIn   ( strCount );
  Vector      strain     ( strCount );

  vector<Vector> gradSlip   ( nslips_ );

  // internal force vector:
  // displacement, slip, tau

  Vector         elemForce1 ( dispCount );  // disp
  vector<Vector> elemForce2 ( nslips_   );  // slip
  vector<Vector> elemForce3 ( nslips_   );  // tau

  // element stiffness matrices 
  // (9 components )

  // Displacement row

  Matrix          elemMat11  ( dispCount, dispCount ); // disp-disp
  vector<Matrix>  elemMat12  ( nslips_ );
  vector<Matrix>  elemMat13  ( nslips_ );

  // Slip rows

  vector<Matrix>  elemMat21  ( nslips_ );
  vector<Matrix>  elemMat22  ( nslips_ );
  vector<Matrix>  elemMat23  ( nslips_ );

  // Schmid stress (tau) rows

  vector<Matrix>  elemMat31  ( nslips_ );
  vector<Matrix>  elemMat32  ( nslips_ );
  vector<Matrix>  elemMat33  ( nslips_ );

  // Resize inner vectors and matrices
  // i.e., corresponding to each slip and tau

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    // Element slip and tau 

    slip[islip] .resize( nodeCount );
    slip0[islip].resize( nodeCount );
    tau[islip]  .resize( nodeCount );

    // Gradients

    gradSlip[islip].resize( rank_ );

    // Element slip and tau internal forces

    elemForce2[islip].resize( nodeCount );
    elemForce3[islip].resize( nodeCount );

    // Element slip and tau matrices

    elemMat12[islip].resize ( dispCount, nodeCount );
    elemMat13[islip].resize ( dispCount, nodeCount );

    elemMat21[islip].resize ( nodeCount, dispCount );
    elemMat22[islip].resize ( nodeCount, nodeCount );
    elemMat23[islip].resize ( nodeCount, nodeCount );

    elemMat31[islip].resize ( nodeCount, dispCount );
    elemMat32[islip].resize ( nodeCount, nodeCount );
    elemMat33[islip].resize ( nodeCount, nodeCount );
  }
  
  // B matrices
  // displacement and scalar fields (slip,tau)
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, nodeCount);
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   dispDofs   ( dispCount );

  vector<IdxVector> slipDofs ( nslips_ );
  vector<IdxVector> tauDofs  ( nslips_ );

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    // Element slip and tau dofs 

    slipDofs[islip].resize( nodeCount );
    tauDofs[islip] .resize( nodeCount );
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

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Get element slip and tau dofs

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      // Element slip and tau dofs 

      dofs_->getDofIndices ( slipDofs[islip], inodes, slipTypes_[islip] );
      dofs_->getDofIndices ( tauDofs[islip],  inodes, tauTypes_[islip]  );
    }

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the shape function

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current element displacements, slips and taus

    disp   = select ( state,  dispDofs );

    // Get element slip and tau dofs

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      // Element slip and tau

      slip[islip]   = select ( state,  slipDofs[islip] );
      slip0[islip]  = select ( state0, slipDofs[islip] );
      tau[islip]    = select ( state,  tauDofs[islip]  );
    }

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      // Element slip and tau internal forces set to zero

      elemForce2[islip] = 0.0;
      elemForce3[islip] = 0.0;
    }

    if ( mbuilder != nullptr )
    {
      elemMat11   = 0.0;

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        // Element slip and tau matrices set to zero

        elemMat12[islip]   = 0.0;
        elemMat13[islip]   = 0.0;

        elemMat21[islip]   = 0.0;
        elemMat22[islip]   = 0.0;
        elemMat23[islip]   = 0.0;

        elemMat31[islip]   = 0.0;
        elemMat32[islip]   = 0.0;
        elemMat33[islip]   = 0.0;
      }    
    }

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {

      // Get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );

      // Compute the displacement B-matrix for this integration point.

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // Get B-matrix associated with phase-field dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the strain for this integration point.

      matmul ( strain,  bd, disp  );

      // Store the regular strain components.

      strain_(slice(BEGIN,strCount),ipoint) = strain;

      // Compute integration point slips, taus, grad slips

      const Vector Nip    = N ( ALL, ip );

      Vector ipSlip  ( nslips_ );
      Vector ipSlip0 ( nslips_ ); 
      Vector ipTau   ( nslips_ );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        // Get integration point slip and tau

        ipSlip[islip]   = dot ( Nip, slip[islip]  );
        ipSlip0[islip]  = dot ( Nip, slip0[islip] );
        ipTau[islip]    = dot ( Nip, tau[islip]   );

        // Get gradient of the slip

        matmul ( gradSlip[islip], be, slip[islip] );
      }  

      // Compute the elastic stress

      material_->update ( stress, stiff, strain, 0  );

      // Compute the inelastic stress

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        stressIn = - ipSlip[islip] * Dsm_[islip];
      }

      // Compute PhiStar (linear)

      Vector ipdPhi ( nslips_ );
      Vector ipPhi  ( nslips_ );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        ipPhi[islip]  =   sign( ipTau[islip] ) / tstar_[islip]
                        * ::pow( fabs( ipTau[islip] ) / 
                                  tauY_[islip], n_[islip] );
        
        ipdPhi[islip] =   sign( ipTau[islip] ) / tstar_[islip]
                        * n_[islip] / ipTau[islip]
                        * ::pow( fabs( ipTau[islip] ) / 
                                  tauY_[islip], n_[islip] );
        
        if ( jem::Float::isNaN( ipdPhi[islip] ) )
        {
          ipdPhi[islip] = 0.0;
        }                        
      }

      // Compute weight of ip

      wip         = ipWeights[ip];

      // -----------------------------------------
      // Compute the internal force components
      // -----------------------------------------

      elemForce1 += wip * mc1.matmul( bdt, Vector(stress - stressIn) );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        elemForce2[islip] += wip * ( Nip * (- dot(Dsm_[islip],strain)
                                            + Esm_(islip,islip) * ipSlip[islip]
                                            + ipTau[islip] ) + 
                                    0.1 * Esm_(islip,islip) * mc1.matmul (bet, gradSlip[islip]) );

        // Note: Esm_ cross hardening is not considered here!

        elemForce3[islip] += wip * ( Nip / dtime_ * ( ipSlip[islip] - ipSlip0[islip] ) - ipPhi[islip] );
      }

      // -----------------------------------------
      // Compute the stiffness matrix components
      // -----------------------------------------

      // Displacement row

      elemMat11  += wip * mc3.matmul ( bdt, stiff, bd );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        // Displacement row

        elemMat12[islip] -= wip * mc2.matmul ( bdt, matmul (Dsm_[islip], Nip ) );
        //elemMat13 = 0.0;  13 component is zero

        // Slip row

        elemMat21[islip] -= wip * mc2.matmul ( matmul (Nip, Dsm_[islip] ), bd );
        elemMat22[islip] += wip * (  matmul( Nip, Nip ) * Esm_(islip,islip)
                               + mc2.matmul ( bet, be ) * Esm_(islip,islip) * 0.1 );
        elemMat23[islip] += wip * matmul( Nip, Nip );

        // Schmid stress (tau) row

        elemMat32[islip] += wip * matmul( Nip, Nip ) * ( 1.0/dtime_ );
        elemMat33[islip] -= wip * matmul( Nip, Nip ) * ipdPhi[islip];
      }

    }  // End of loop on integration points

    // Assembly ...

    if ( mbuilder != nullptr )
    {
      mbuilder -> addBlock ( dispDofs, dispDofs, elemMat11 );

      for ( int islip = 0; islip < nslips_; islip++ )
      {
        mbuilder -> addBlock ( dispDofs, slipDofs[islip], elemMat12[islip] );

        mbuilder -> addBlock ( slipDofs[islip], dispDofs,        elemMat21[islip] );
        mbuilder -> addBlock ( slipDofs[islip], slipDofs[islip], elemMat22[islip] );
        mbuilder -> addBlock ( slipDofs[islip], tauDofs[islip],  elemMat23[islip] );

        mbuilder -> addBlock ( tauDofs[islip], slipDofs[islip],  elemMat32[islip] );
        mbuilder -> addBlock ( tauDofs[islip], tauDofs[islip],   elemMat33[islip] );
      }      
    }

    select ( force, dispDofs ) += elemForce1;

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      select ( force, slipDofs[islip] ) += elemForce2[islip];
      select ( force, tauDofs[islip]  ) += elemForce3[islip];
    }
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void CrystalViscoPlasticityModel::getMatrix2_

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

void CrystalViscoPlasticityModel::getArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  CrystalViscoPlasticityModel::initializeIPMPMap_ ( )

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


bool CrystalViscoPlasticityModel::getTable_

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
    getStress_ ( *table, weights);

    return true;
  }

  // Nodal strains
  if ( name == "strain" && 
       table->getRowItems() == nodes_.getData() )
  {
    getStrain_ ( *table, weights);

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


void CrystalViscoPlasticityModel::getStress_

  ( XTable&        table,
    const Vector&  weights )

{
}


//-----------------------------------------------------------------------
//   getStrain_
//-----------------------------------------------------------------------


void CrystalViscoPlasticityModel::getStrain_

  ( XTable&        table,
    const Vector&  weights )

{
}


//-----------------------------------------------------------------------
//   getOutputData_
//-----------------------------------------------------------------------

/** all data at Gauss points are computed, however, only those requested
 *  in the vtk block in the .pro file will be written to vtu files.
 *  Example: vtk.data = "stress_yy | stress_xx" for xx and yy stress 
 *  components.
 */

void CrystalViscoPlasticityModel::getOutputData_

  ( Ref<XTable>    table,
    const Vector&  weights,
    const String&  contents,
    const Vector&  state )
          
{
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void CrystalViscoPlasticityModel::checkCommit_

  ( const Properties&  params )

{
  // System::info() << myName_ << " : check commit ... do nothing!\n";
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newCrystalViscoPlasticityModel
//-----------------------------------------------------------------------


static Ref<Model>     newCrystalViscoPlasticityModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<CrystalViscoPlasticityModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareCrystalViscoPlasticityModel
//-----------------------------------------------------------------------


void declareCrystalViscoPlasticityModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "CrystalViscoPlasticity", & newCrystalViscoPlasticityModel );
}
