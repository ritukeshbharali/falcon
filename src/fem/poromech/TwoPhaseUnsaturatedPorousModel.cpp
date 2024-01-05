
/** @file TwoPhaseUnsaturatedPorousModel.cpp
 *  @brief Two-phase unsaturated porous media model.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 22 April 2022
 * 
 *  NOTE: Stiffness matrix is unsymmetric.
 *
 *  Updates (when, what and who)
 *     - [19 May 2022] Corrected the expression for the
 *       storage term, Sto (RB) [BUG FIX!]
 *     - [19 October 2022] Updated to a mass conserving
 *       scheme (RB)
 *     - [12 December 2022] Retention models update
 *       similar to material models (RB)
 *     - [27 December 2023] getMatrix_ returns the internal
 *       force if mbuilder = nullptr. Eliminates duplicate
 *       code. (RB)  
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/geom/Geometries.h>

/* Include other headers */

#include "FalconPoroMechModels.h"
#include "TwoPhaseUnsaturatedPorousModel.h"
#include "util/TbFiller.h"
#include "util/XNames.h"
#include "util/Constants.h"


//=======================================================================
//   class TwoPhaseUnsaturatedPorousModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  TwoPhaseUnsaturatedPorousModel::DISP_NAMES[3]        = { "dx", "dy", "dz" };
const char*  TwoPhaseUnsaturatedPorousModel::SHAPE_PROP           = "shape";
const char*  TwoPhaseUnsaturatedPorousModel::MATERIAL_PROP        = "material";
const char*  TwoPhaseUnsaturatedPorousModel::RETENTION_PROP       = "retention";

const char*  TwoPhaseUnsaturatedPorousModel::INTRIN_PERM_PROP     = "intrinPerm";
const char*  TwoPhaseUnsaturatedPorousModel::FLUID_VISC_PROP      = "fluidVisc";
const char*  TwoPhaseUnsaturatedPorousModel::SOLID_STIFF_PROP     = "solidStiff";
const char*  TwoPhaseUnsaturatedPorousModel::FLUID_STIFF_PROP     = "fluidStiff";
const char*  TwoPhaseUnsaturatedPorousModel::POROSITY_PROP        = "porosity";
const char*  TwoPhaseUnsaturatedPorousModel::BIOT_COEFF_PROP      = "biotCoeff";
const char*  TwoPhaseUnsaturatedPorousModel::DTIME_PROP           = "dtime";

const char*  TwoPhaseUnsaturatedPorousModel::SOLID_DENSITY_PROP   = "rhoSolid";
const char*  TwoPhaseUnsaturatedPorousModel::FLUID_DENSITY_PROP   = "rhoFluid";
const char*  TwoPhaseUnsaturatedPorousModel::GRAVITY_PROP         = "gravity";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


TwoPhaseUnsaturatedPorousModel::TwoPhaseUnsaturatedPorousModel

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


  // Get the DOF space, add displacement and pressure DOFs

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ + 1 );

  // Make a shallow copy of dofTypes_ for displacement.

  dispTypes_.ref ( dofTypes_[slice(BEGIN,rank_)] );

  for ( int i = 0; i < rank_; i++ )
  {
    dispTypes_[i] = dofs_->addType ( DISP_NAMES[i] );
  }

  // Make a shallow copy of dofTypes_ for pressure.

  presTypes_.ref ( dofTypes_[slice(rank_,END)] );

  presTypes_[0] = dofs_->addType ( "dp" );
    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  // Allocate memory for the strains and set them to zero

  strain_.resize ( STRAIN_COUNTS[rank_], ipCount );
  strain_ = 0.0;

  // Create a material model object.

  material_ = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  material_-> allocPoints  ( ipCount );

  initializeIPMPMap_ ();

  // Create a retention model object.

  retention_ = newRetentionMaterial ( RETENTION_PROP, myConf, myProps, globdat );

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  IdxVector   ielems     = egroup_.getIndices ();
  const int   ielemCount = ielems.size        ();

  // Set all elements to active.

  isActive_.resize ( ielemCount );
  isActive_ = 1;


  // Configure poromechanical parameters

  kappa_       = 1.0e-12;
  myProps.find ( kappa_, INTRIN_PERM_PROP );
  myConf. set  ( INTRIN_PERM_PROP, kappa_ );

  mu_          = 8.9e-4;
  myProps.find ( mu_, FLUID_VISC_PROP );
  myConf. set  ( FLUID_VISC_PROP, mu_ );

  Ks_          = 1.0e+10;
  myProps.find ( Ks_, SOLID_STIFF_PROP );
  myConf. set  ( SOLID_STIFF_PROP, Ks_ );

  Kf_          = 2.0e+09;
  myProps.find ( Kf_, FLUID_STIFF_PROP );
  myConf. set  ( FLUID_STIFF_PROP, Kf_ );

  phi_         = 0.375;
  myProps.find ( phi_, POROSITY_PROP );
  myConf. set  ( POROSITY_PROP, phi_ );

  alpha_      = 1.0;
  myProps.find ( alpha_, BIOT_COEFF_PROP );
  myConf. set  ( BIOT_COEFF_PROP, alpha_ );

  dtime_       = 1.0;
  myProps.find ( dtime_, DTIME_PROP );
  myConf. set  ( DTIME_PROP, dtime_ );

  rhoS_        = 2000.0;
  myProps.find ( rhoS_, SOLID_DENSITY_PROP );
  myConf. set  ( SOLID_DENSITY_PROP, rhoS_ );

  rhoF_        = 1000.0;
  myProps.find ( rhoF_, FLUID_DENSITY_PROP );
  myConf. set  ( FLUID_DENSITY_PROP, rhoF_ );

  // Compute derived quantities

  gravity_     = true;
  myProps.find ( gravity_, GRAVITY_PROP );
  myConf. set  ( GRAVITY_PROP, gravity_ );

  if (gravity_)
  {
    g_    = 9.81;
  }
  else
  {
    g_    = 0.0;
  }
  
  Keff_ = kappa_ / mu_ ;

  gVec_.resize ( rank_ );
  gVec_ = 0.0;
  gVec_ [rank_-1] = -g_;


  // Set up Voigt divergence operator based of mesh rank

  voigtDiv_.resize ( STRAIN_COUNTS[rank_]);
  voigtDiv_ = 0.0;

  if ( rank_ == 1 )
  {
    voigtDiv_[0] = 1.0;      // xx-component
  }
  else if ( rank_ == 2 )
  {
    voigtDiv_[0] = 1.0;      // xx-component
    voigtDiv_[1] = 1.0;      // yy-component
  }
  else if ( rank_ == 3 )
  {
    voigtDiv_[0] = 1.0;      // xx-component
    voigtDiv_[1] = 1.0;      // yy-component
    voigtDiv_[2] = 1.0;      // zz-component
  }
  else
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  }
 
}


TwoPhaseUnsaturatedPorousModel::~TwoPhaseUnsaturatedPorousModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );
  Properties  retProps = myProps.findProps ( RETENTION_PROP );

  material_-> configure ( matProps, globdat );
  retention_->configure ( retProps, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP  );
  Properties  retConf = myConf.makeProps ( RETENTION_PROP );

  material_-> getConfig ( matConf, globdat );
  retention_->getConfig ( retConf, globdat );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool TwoPhaseUnsaturatedPorousModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // compute the internal force vector

  if ( action == Actions::GET_INT_VECTOR )
  {
    Vector  state;
    Vector  state0;
    Vector  force;

    // Get the current and old step displacements.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( force, ActionParams::INT_VECTOR );

    getMatrix_ ( nullptr, force, state, state0 );

    return true;
  }

  if ( action == Actions::GET_MATRIX0 )
  {

    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  state0;
    Vector  force;

    // Get the current and old step displacements.

    StateVector::get    ( state, dofs_, globdat );
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

  /** During initialization of the model, the pressure dofs of the
   *  mid nodes of all domain elements are linearly constrained to
   *  the corner nodes of the respective edge. This yields a Taylor
   *  Hood element.
  */ 

  if ( action == Actions::INIT )
  {
    // Get the constraints associated with the DOF space.

    Ref<Constraints>  cons = Constraints::get ( dofs_, globdat );

    setConstraints_ ( *cons );

    return true;
  }

  /** COMMIT: Requests an action when a computation step has converged.
   *  At this point, the material internal state variables are updated.
   */ 

  if ( action == Actions::COMMIT )
  {
    material_->commit ();

    return true;
  }

  /** CHECK_COMMIT: Can used to discard the current step
  */ 

  if ( action == Actions::CHECK_COMMIT )
  {
    // checkCommit_ ( params );
    return true;
  }

  /** GET_TABLE: Requests a post-processing action, such as computing
   *  tabular data (stress, strains, etc.) to print VTK files.
   */ 
  
  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  /** SET_STEP_SIZE: Implements changes to the step-size, thrown by a
   *  solver module (e.g., AdaptiveStepping, ReduceStepping)
   */ 

  if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::getMatrix_

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

  // number of displacement and pressure dofs

  const int   dispCount  = nodeCount * rank_;  
  const int   presCount  = nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement and pressure
  
  Vector      disp       ( dispCount );
  Vector      pres       ( presCount );

  // old step element vector state:
  // displacement and pressure
  
  Vector      disp0      ( dispCount );
  Vector      pres0      ( presCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );

  // internal force vector:
  // displacement and pressure

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( presCount  );

  // element stiffness matrices (four components)

  Matrix      elemMat1   ( dispCount, dispCount ); // disp-disp
  Matrix      elemMat2   ( dispCount, presCount ); // disp-pres
  Matrix      elemMat3   ( presCount, dispCount ); // pres-disp
  Matrix      elemMat4   ( presCount, presCount ); // pres-pres
  
  // B matrices
  // displacement and pressure
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, presCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   dispDofs   ( dispCount );
  IdxVector   presDofs   ( presCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Initialize saturation, relative permeability and their derivatives

  double Sf      = 0.0;
  double krw     = 0.0;
  double dSfdp   = 0.0;
  double dkrwdp  = 0.0;

  double Sf0      = 0.0;
  double krw0     = 0.0;
  double dSfdp0   = 0.0;
  double dkrwdp0  = 0.0;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {

    if ( ! isActive_[ie] )
    {      
      continue;
    }

    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );
    dofs_->getDofIndices ( presDofs, inodes, presTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current nodal displacements and pressures

    disp = select ( state, dispDofs );
    pres = select ( state, presDofs );

    // Get old step nodal displacements and pressures

    disp0 = select ( state0, dispDofs );
    pres0 = select ( state0, presDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0;

    if ( mbuilder != nullptr )
    {
      elemMat1   = 0.0;
      elemMat2   = 0.0;
      elemMat4   = 0.0;
    }

    elemMat3   = 0.0;

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {

      // get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );

      // Compute the B-matrix for this integration point.
      // it is the B-matrix of displacement dofs

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // get B-matrix associated with pres dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the strain for this integration point.

      matmul ( strain,  bd, disp  );

      // Store the regular strain components.

      strain_(ALL,ipoint) = strain;

      // Compute the integration point pressures

      Vector Nip          = N( ALL,ip );

      double p            = dot ( Nip, pres  );
      double p0           = dot ( Nip, pres0 );

      // Update the material model (get stress and stiffness).

      material_->update ( stress, stiff, strain_(ALL,ipoint), ipoint );

      // Update the retention model (get saturation, rel. permeability, derivatives).

      retention_->update( Sf,  dSfdp,  krw,  dkrwdp,  p  );
      retention_->update( Sf0, dSfdp0, krw0, dkrwdp0, p0 );


      // Compute the storage coefficient and its derivative

      double Sto1    = ( ( alpha_ - phi_ ) / Ks_ * Sf * Sf ) + ( phi_ * Sf / Kf_ );
      double Sto2    = ( ( alpha_ - phi_ ) / Ks_ * Sf * p ) + phi_;

      double dSto1dp = ( Sf / Ks_ * ( alpha_ - phi_ ) * 2.0 * dSfdp ) + ( phi_ * dSfdp / Kf_ );
      double dSto2dp = ( ( alpha_ - phi_ ) / Ks_ * Sf ) + ( ( alpha_ - phi_ ) / Ks_ * dSfdp * p);

      // compute weight of the ip

      wip         = ipWeights[ip];

      // compute stiffness matrix components

      if ( mbuilder != nullptr )
      {
        elemMat1   += wip * mc3.matmul ( bdt, stiff, bd );
        elemMat2   -= wip * alpha_ * Sf * mc2.matmul ( bdt, matmul (voigtDiv_, Nip ) );
        elemMat2   -= wip * alpha_ * dSfdp * p * mc2.matmul ( bdt, matmul (voigtDiv_, Nip ) );
        elemMat4   += wip * ( ( Sto1 + Sto2 * dSfdp + dSto1dp * ( p-p0 ) 
                               + dSto2dp * ( Sf-Sf0 ) ) * matmul ( Nip, Nip ) 
                               + dtime_ * krw * Keff_ * mc2.matmul ( bet, be ) );
      }

      elemMat3   += wip * alpha_ * Sf * mc2.matmul ( matmul ( Nip, voigtDiv_ ), bd );
           
      // compute internal forces

      elemForce1 +=  wip * ( mc1.matmul ( bdt, stress ) );
      elemForce1 -=  wip * alpha_ * Sf * p * ( mc1.matmul ( bdt, voigtDiv_ ) );
      elemForce2 +=  wip * ( (Sto1 * (p-p0) + Sto2 * (Sf - Sf0) ) * Nip 
                     + dtime_ * krw * Keff_ * mc1.matmul( mc2.matmul ( bet, be ), pres ) );
      elemForce2 -=  wip * dtime_ * krw * Keff_ * rhoF_ * mc1.matmul( bet, gVec_ );

      // compute gravity body load ( -y or -z direction for 2D/3D)

      elemForce1[slice(rank_-1,END,rank_)] -= Nip * wip * ( phi_ * Sf * rhoF_ + ( 1.0 - phi_ ) * rhoS_ ) * gVec_[rank_-1];

    }  // end of loop on integration points

    elemForce2   +=  mc1.matmul ( elemMat3, disp ) - mc1.matmul ( elemMat3, disp0 );

    // Assembly ...

    if ( mbuilder != nullptr )
    {
      mbuilder -> addBlock ( dispDofs, dispDofs, elemMat1 );
      mbuilder -> addBlock ( dispDofs, presDofs, elemMat2 );
      mbuilder -> addBlock ( presDofs, dispDofs, elemMat3 );
      mbuilder -> addBlock ( presDofs, presDofs, elemMat4 );
    }

    select ( force, dispDofs ) += elemForce1;
    select ( force, presDofs ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void TwoPhaseUnsaturatedPorousModel::getMatrix2_

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
    nodes_.getSomeCoords ( coords, inodes   );
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

    // Assemble mass matrix

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}

//-----------------------------------------------------------------------
//   setConstraints_
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::setConstraints_ ( Constraints& cons )
{
  using jive::util::Printer;
  using jive::geom::Geometries;

  IdxVector     ielems     = egroup_.getIndices ();

  const int     ielemCount = ielems   .size     ();
  const int     nodeCount  = shape_->nodeCount  ();

  const String  shapeGeom  = shape_->getGeometry ();

  IdxVector     inodes     ( nodeCount );
  IdxVector     idofs      ( nodeCount );
  IdxVector     jdofs      ( 2 );
  Vector        coeffs     ( 2 );

  IntMatrix     cdofs;


  // only do this constraint for high order elements
  // not for linear elements

  if      ( shapeGeom == Geometries::LINE &&
      nodeCount == 2 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::TRIANGLE &&
      nodeCount == 3 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::SQUARE &&
      nodeCount == 4 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::TETRAHEDRON &&
      nodeCount == 4 )
  {
    return;
  }

  System::info() << "Two Phase Unsaturated Porous Model ...\n"
     << "   : Initializing constraints for Taylor-Hood element ...\n";

  // Determine which pressure DOFS are to be constrained. The
  // first row of the cdofs matrix contains the DOF indices to be
  // constrained, and the other two rows contain the DOF indices of
  // the two master nodes. Note that the DOF indices are local with
  // respect to an element. They will be translated to global DOF
  // indices below.

  if      ( shapeGeom == Geometries::LINE &&
      nodeCount == 3 )
  {
    cdofs.resize ( 3, 1 );

    // Constrain the mid node.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;
  }
  else if ( shapeGeom == Geometries::TRIANGLE &&
      nodeCount == 6 )
  {
    cdofs.resize ( 3, 3 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 0;
  }
  else if ( shapeGeom == Geometries::SQUARE &&
      nodeCount == 8 )
  {
    cdofs.resize ( 3, 4 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 6;

    cdofs(0,3) = 7;
    cdofs(1,3) = 6;
    cdofs(2,3) = 0;
  }
  else
  {
    throw IllegalInputException (
      getContext (),
      String::format (
  "unsupported shape geometry: %S%d",
  &shapeGeom,
  nodeCount
      )
    );
  } 

  // Determine the constraint coefficients.

  coeffs = 0.5;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element DOFs.

    elems_.getElemNodes  ( inodes, ielem );
    dofs_->getDofIndices ( idofs,  inodes, presTypes_ );

    // Add constraints.

    for ( int j = 0; j < cdofs.size(1); j++ )
    {
      // Get the global index of the DOF to be constrained.

      int  idof = idofs[cdofs(0,j)];

      // Check whether this DOF has already been constrained.

      if ( cons.isSlaveDof( idof ) )
      {
        continue;
      }

      // Get the global indices of the two master DOFs.

      jdofs[0] = idofs[cdofs(1,j)];
      jdofs[1] = idofs[cdofs(2,j)];

      cons.addConstraint ( idof, jdofs, coeffs );
    }
  }
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  TwoPhaseUnsaturatedPorousModel::initializeIPMPMap_ ( )

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


bool TwoPhaseUnsaturatedPorousModel::getTable_

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
    Vector disp;

    StateVector::get (disp, dofs_, globdat);

    getStress_ ( *table, weights);

    return true;
  }

  // Nodal strains
  if ( name == "strain" && 
       table->getRowItems() == nodes_.getData() )
  {
    Vector disp;

    StateVector::get (disp, dofs_, globdat);

    getStrain_ ( *table, weights);

    return true;
  }

  // 
  if ( name == "xoutTable" )
  {
    Vector  disp;
    String  contents;

    StateVector::get ( disp,     dofs_, globdat  );
    params.      get ( contents, "contentString" );

    
    getOutputData_ ( table, weights, contents, disp );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::getStress_

  ( XTable&        table,
    const Vector&  weights )

{
  Matrix     sfuncs     = shape_->getShapeFunctions ();

  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];

  Matrix     stiff      ( strCount, strCount  );
  Matrix     ndStress   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  Vector     stressIp   ( strCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  MChain1    mc1;

  int        ipoint, igpoint = 0;

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
    int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount; ip++, igpoint++ )
    {
      // Extrapolate the integration point stresses to the nodes using
      // the transposed shape functions.

      ipoint     = ipMpMap_ (ielem,ip);

      material_-> update ( stressIp, stiff, strain_(ALL,igpoint) , ipoint );

      ndStress  += matmul ( sfuncs(ALL,ip), stressIp );
      ndWeights += sfuncs(ALL,ip);
    }

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols, ndStress );
  }
}


//-----------------------------------------------------------------------
//   getStrain_
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::getStrain_

  ( XTable&        table,
    const Vector&  weights )

{
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount = ielems.size         ();
  const int  nodeCount  = shape_->nodeCount   ();
  const int  ipCount    = shape_->ipointCount ();
  const int  strCount   = STRAIN_COUNTS[rank_];

  Matrix     sfuncs     = shape_->getShapeFunctions ();

  Matrix     ndStrain   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  int        ipoint;


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

  ipoint = 0;

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element nodes.

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndStrain  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
    {
      // Extrapolate the integration point strains to the nodes using
      // the transposed shape functions.

      ndStrain  += matmul ( sfuncs(ALL,ip), strain_(ALL,ipoint) );
      ndWeights += sfuncs(ALL,ip);
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
//   getHistory_
//-----------------------------------------------------------------------


void TwoPhaseUnsaturatedPorousModel::getHistory_

  ( XTable&          table,
    const Vector&    weights )
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

void TwoPhaseUnsaturatedPorousModel::getOutputData_

  ( Ref<XTable>    table,
    const Vector&  weights,
    const String&  contents,
    const Vector&  state )
          
{
  MChain2    mc2;

  // table filler related stuff

  TbFiller   tbFiller   ( rank_ );

  //StringVector hisNames =  materials_[0]->getHistoryNames ();

  Slice      iistrain   = tbFiller.announce ( "strain.tensor" ); // iistrain  = [0 1 2 3]
  Slice      iistress   = tbFiller.announce ( "stress.tensor" ); // iistrain  = [0 1 2 3]

  Vector     ipValues   ( tbFiller.typeCount() ); // typeCount() = # of types = 8 in 2D

  Vector     strain     ( ipValues[iistrain]   );
  Vector     stress     ( ipValues[iistress]   );

  // Let TbFiller find out which columns of ndValues to write to 
  // which columns of the table (based on filter in input file)

  IdxVector  i2table;
  IdxVector  jcols;

  tbFiller . setFilter   ( contents ); 
  tbFiller . prepareTable( i2table, jcols, table ); 
  
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];
  const int  dispCount   = nodeCount * rank_;
   
  Matrix     ndValuesOut  ( nodeCount, i2table.size() );
  Matrix     ndValuesOut1 ( nodeCount, i2table.size() );
  Vector     ipValuesOut  ( i2table.size() );

  Matrix      stiff      ( strCount, strCount  );
  Matrix      bd         ( strCount, dispCount );
  Vector      ndWeights  ( nodeCount           ); 
  IdxVector   inodes     ( nodeCount           );


  Properties  params;

    idx_t       ipoint, igpoint = 0;

  Matrix      sfuncs     = shape_->getShapeFunctions ();

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    ndValuesOut = 0.;
    ndWeights   = 0.;

    // Loop on integration points 
    
    for ( idx_t ip = 0; ip < ipCount; ip++, igpoint++ )
    {
      ipoint      = ipMpMap_( ielem, ip );
      
 
      strain      = strain_(slice(BEGIN,strCount),igpoint);
      material_-> update ( stress, stiff, strain , ipoint );

       // apply the filter now, only cols specified by i2table are 
       // written to the table 
      ipValuesOut  = ipValues[i2table]; 

      matmul (ndValuesOut1, sfuncs(ALL,ip), ipValuesOut );


      ndValuesOut += ndValuesOut1;
      ndWeights   += sfuncs(ALL,ip);

    }  // end of loop on integration points
    
    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table->addBlock ( inodes, jcols, ndValuesOut );
  }
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void TwoPhaseUnsaturatedPorousModel::checkCommit_

  ( const Properties&  params )

{
  System::info() << myName_ << " : check commit ... do nothing!\n";
 
}

//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void TwoPhaseUnsaturatedPorousModel::setStepSize_

  ( const Properties&  params )

{

  System::out() << "Enters setStepSize_ \n";

  double       dt;

  params.get ( dt,  XProps::STEP_SIZE   );

  System::out() << "Setting step size to " << dt << " ("
    << dt / dtime_ * 100 << "\% of previous step size)" << "\n";

  dtime_ = dt;
  
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newTwoPhaseUnsaturatedPorousModel
//-----------------------------------------------------------------------


static Ref<Model>     newTwoPhaseUnsaturatedPorousModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<TwoPhaseUnsaturatedPorousModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareTwoPhaseUnsaturatedPorousModel
//-----------------------------------------------------------------------


void declareTwoPhaseUnsaturatedPorousModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "TwoPhaseUnsaturatedPorous", & newTwoPhaseUnsaturatedPorousModel );
}
