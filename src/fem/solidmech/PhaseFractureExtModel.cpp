
/** @file PhaseFractureExtModel.cpp
 *  @brief Implements a extrapolation-based convexified phase-field 
 *         fracture model.
 *  
 *  This class implements the unified phase-field fracture
 *  model (see DOI: 10.1016/j.jmps.2017.03.015). Fracture
 *  irreversibility is enforced using the penalization 
 *  approach (see DOI: 10.1016/j.cma.2019.05.038). The
 *  fracture driving and resisting energies are obtained
 *  using Amor split (see DOI: 10.1016/j.jmps.2009.04.011).
 *  Convexification is carried out based on extrapolation
 *  of the phase-field for the momentum balance equation
 *  (see DOI: 10.1016/j.cma.2015.03.009).
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 30 March 2022
 * 
 *  @note This model must be run with ReduceStepping
 *        Module.
 *  
 *  Updates (when, what and who)
 *     - [14 June 2022], added getIntForce_ and removed
 *       getHistory_. getIntForce_ is required for line
 *       search operations. (RB)
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>

/* Include other headers */

#include "FalconSolidMechModels.h"
#include "PhaseFractureExtModel.h"
#include "util/TbFiller.h"
#include "util/XNames.h"

#include "materials/Material.h"
#include "materials/PhaseFractureMaterial.h"


//=======================================================================
//   class PhaseFractureExtModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  PhaseFractureExtModel::DISP_NAMES[3]         = { "dx", "dy", "dz" };
const char*  PhaseFractureExtModel::SHAPE_PROP            = "shape";
const char*  PhaseFractureExtModel::MATERIAL_PROP         = "material";
const char*  PhaseFractureExtModel::RHO_PROP              = "rho";

const char*  PhaseFractureExtModel::FRACTURE_TYPE_PROP    = "fractureType";
const char*  PhaseFractureExtModel::GRIFFITH_ENERGY_PROP  = "griffithEnergy";
const char*  PhaseFractureExtModel::LENGTH_SCALE_PROP     = "lengthScale";
const char*  PhaseFractureExtModel::TENSILE_STRENGTH_PROP = "tensileStrength";
const char*  PhaseFractureExtModel::PENALTY_PROP          = "penalty";

const char*  PhaseFractureExtModel::BRITTLE_AT1           = "at1";
const char*  PhaseFractureExtModel::BRITTLE_AT2           = "at2";
const char*  PhaseFractureExtModel::LINEAR_CZM            = "czmLinear";
const char*  PhaseFractureExtModel::EXPONENTIAL_CZM       = "czmExponential";
const char*  PhaseFractureExtModel::CORNELISSEN_CZM       = "czmCornelissen";
const char*  PhaseFractureExtModel::HYPERBOLIC_CZM        = "czmHyperbolic";

vector<String> PhaseFractureExtModel::fractureModels (6);


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


PhaseFractureExtModel::PhaseFractureExtModel

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

  // Make a shallow copy of the first part of the dofTypes_ array.

  dispTypes_.ref ( dofTypes_[slice(BEGIN,rank_)] );

  for ( int i = 0; i < rank_; i++ )
  {
    dispTypes_[i] = dofs_->addType ( DISP_NAMES[i] );
  }

  // Make a shallow copy of the last part of the dofTypes_ array.

  phiTypes_.ref ( dofTypes_[slice(rank_,END)] );

  phiTypes_[0] = dofs_->addType ( "phi" );
    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  // Allocate memory for the strains and set them to zero

  strain_.resize ( STRAIN_COUNTS[rank_] + 1, ipCount );
  strain_ = 0.0;

  // Create a material model object. Ensure it is a phase-field material 
  // with dynamicCast

  material_      = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  phaseMaterial_ = dynamicCast<PhaseFractureMaterial> ( material_ );
  
  material_->allocPoints  ( ipCount );

  initializeIPMPMap_ ();

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  IdxVector   ielems     = egroup_.getIndices ();
  const int   ielemCount = ielems.size         ();

  isActive_.resize ( ielemCount );
  isActive_ = 1;

  /*
   *  Setup the fracture models
   */

  fractureModels[0]     = BRITTLE_AT1;
  fractureModels[1]     = BRITTLE_AT2;
  fractureModels[2]     = LINEAR_CZM;
  fractureModels[3]     = EXPONENTIAL_CZM;
  fractureModels[4]     = CORNELISSEN_CZM;
  fractureModels[5]     = HYPERBOLIC_CZM;

  // Get fracture type

  myProps.get ( fracType_, FRACTURE_TYPE_PROP  );

  if ( std::find ( fractureModels.begin (),
                   fractureModels.end   (),
                   fracType_ ) == fractureModels.end () )
  {
    throw Error (
      JEM_FUNC,
      String("unexpected definition of fracture type!!!\n") +
      String("Supported types include: \n") +
      BRITTLE_AT1   + String(", ") + BRITTLE_AT2 + String(", ") +
      LINEAR_CZM + String(", ") + EXPONENTIAL_CZM + String(", ") +
      CORNELISSEN_CZM + String(", ") + HYPERBOLIC_CZM
    );
  }

  myConf. set  ( FRACTURE_TYPE_PROP, fracType_ );

  // Configure phase-field fracture parameters

  gc_          = 2.7;
  myProps.find ( gc_, GRIFFITH_ENERGY_PROP );
  myConf. set  ( GRIFFITH_ENERGY_PROP, gc_ );

  l0_          = 0.015;
  myProps.find ( l0_, LENGTH_SCALE_PROP );
  myConf. set  ( LENGTH_SCALE_PROP, l0_ );

  sigmaT_      = 0.0;
  myProps.find ( sigmaT_, TENSILE_STRENGTH_PROP );
  myConf. set  ( TENSILE_STRENGTH_PROP, sigmaT_ );

  beta_         = 6500.0;
  myProps.find ( beta_, PENALTY_PROP );
  myConf. set  ( PENALTY_PROP, beta_ );

  dt_          = 0.0;
  dt0_         = 0.0;
}


PhaseFractureExtModel::~PhaseFractureExtModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void PhaseFractureExtModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  material_->configure ( matProps, globdat );

  // Configure phase-field fracture parameters

  if ( fracType_ == BRITTLE_AT1 )
  {
    cw_     = 8.0/3.0;
    eta_    = 1.0;
    p_      = 2.0;
    a1_     = 2.0;
    a2_     = -0.5;
    a3_     = 0.0;
    Psi0_   = (3. * gc_)/ (16. * l0_);
  }
  else if ( fracType_ == BRITTLE_AT2 )
  {
    cw_    = 2.0;
    eta_   = 0.0;
    p_     = 2.0;
    a1_    = 2.0;
    a2_    = -0.5;
    a3_    = 0.0;
    Psi0_  = 0.0;
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void PhaseFractureExtModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP );

  material_->getConfig ( matConf, globdat );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool PhaseFractureExtModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // Compute the internal force vector

  if ( action == Actions::GET_INT_VECTOR )
  {

    Vector  state;
    Vector  state0;
    Vector  state00;
    Vector  force; 

    // Get the current and old step displacements.

    StateVector::get       ( state,   dofs_, globdat );
    StateVector::getOld    ( state0,  dofs_, globdat );
    StateVector::getOldOld ( state00, dofs_, globdat );
    
    // Get the internal force vector.

    params.get ( force,    ActionParams::INT_VECTOR );

    getIntForce_ ( force, state, state0, state00 );

    return true;
  }

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  state0;
    Vector  state00;
    Vector  force; 

    // Get the current and old step displacements.

    StateVector::get       ( state,   dofs_, globdat );
    StateVector::getOld    ( state0,  dofs_, globdat );
    StateVector::getOldOld ( state00, dofs_, globdat );
    
    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, state, state0, state00 );

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
    // material_->commit (); // not required for linear elastic material

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


  if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getIntForce_
//-----------------------------------------------------------------------


void PhaseFractureExtModel::getIntForce_

  ( const Vector&   force,
    const Vector&   state,
    const Vector&   state0,
    const Vector&   state00 )

{
  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of displacement and phase-field dofs

  const int   dispCount  = nodeCount * rank_;  
  const int   phiCount   = nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement and phase-field
  
  Vector      disp       ( dispCount );
  Vector      phi        ( phiCount );

  // old step element vector state:
  // phase-field
  
  Vector      phi0       ( phiCount );

  // old old step element vector state:
  // phase-field
  
  Vector      phi00       ( phiCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );
  Vector      stressP    ( strCount );

  // internal force vector:
  // displacement and phase-field

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( phiCount  );
  
  // B matrices
  // displacement and phase-field
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, phiCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   phiDofs    ( phiCount  );
  IdxVector   dispDofs   ( dispCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  double      gphiEx;
  double      dgphi;

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
    dofs_->getDofIndices ( phiDofs, inodes, phiTypes_ );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the shape function

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current element displacements and phase-fields

    phi   = select ( state, phiDofs  );
    disp  = select ( state, dispDofs );

    // Get old step and old old step element phase-field

    phi0  = select ( state0,  phiDofs );
    phi00 = select ( state00, phiDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0; 

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

      // Compute the integration point phase-fields (current, old, old old)

      Vector Nip           = N( ALL,ip );

      double pf            = dot ( Nip, phi   );
      double pf0           = dot ( Nip, phi0  );
      double pf00          = dot ( Nip, phi00 );
      double pfEx          = 0.0;

      // Linearly extrapolate old old phase-field to current step

      if (dt0_ < 1e-20 )
      {
        pfEx = pf0;
      }
      else
      {
        // pfEx              = pf00 + ( dt_ + dt0_ )/dt0_ * ( pf0 - pf00 );
        // corrected Wick book
        pfEx = - pf00 * dt_/dt0_ + pf0 * (dt_+dt0_)/dt0_;
      }
      
      // Ensure extrapolated phase-field is within bounds [0,1]

      if ( pfEx < 0.0 )
      {
        pfEx = 0.0;
      }
      else if ( pfEx > 1.0 )
      {
        pfEx = 1.0;
      }

      // Compute the degradation function with the extrapolated phase-field

      {
        double gphiEx_n    = ::pow( 1- pfEx, p_);
        double gphiEx_d    = gphiEx_n + a1_*pfEx + a1_*a2_*pfEx*pfEx + 
                             a1_*a2_*a3_*pfEx*pfEx*pfEx;
        gphiEx             = gphiEx_n/gphiEx_d + 1.e-7;

      }

      strain_(strCount,ipoint) = gphiEx;

      // Compute the degradation function derivatives with the current phase-field

      {
        double gphi_n    = ::pow( 1- pf, p_);
        double gphi_d    = gphi_n + a1_*pf + a1_*a2_*pf*pf + a1_*a2_*a3_*pf*pf*pf;

        double dgphi_n   = - p_ * ( ::pow( 1- pf, p_- 1) );
        double dgphi_d   = dgphi_n + a1_ + 2.0*a1_*a2_*pf + 3.0*a1_*a2_*a3_*pf*pf;

        dgphi            = ( dgphi_n*gphi_d - gphi_n*dgphi_d ) / ( gphi_d * gphi_d );
      }

      // Compute the locally dissipated fracture energy function derivative

      double dw        = eta_ + 2.0 * ( 1- eta_ ) * pf;

      // Compute stress and material tangent stiffness

      phaseMaterial_->update ( stress, stiff, strain, gphiEx, ipoint );

      // Get fracture driving energy

      double Psi = phaseMaterial_->givePsi();
      Psi        = max(Psi,Psi0_);

      // Compute stiffness matrix components

      wip         = ipWeights[ip];
     
      // Compute internal forces

      elemForce1 +=  wip * ( mc1.matmul ( bdt, stress ) );
      elemForce2 +=  wip * ( Nip * ( gc_/(cw_*l0_) * dw + dgphi * Psi) 
                     + (2.0 * gc_*l0_/cw_) 
                     * mc1.matmul ( mc2. matmul ( bet, be ), phi ) );

      // Penalization (enforce fracture irreversibility) 

      if ( pf0 - pf > 1.e-15 )
      {
        double Dpf = pf0 - pf;
        elemForce2 -= wip * ( gc_/(cw_*l0_) * beta_ * Dpf ) * Nip;
      }   

    }  // End of loop on integration points

    // Assembly ...

    select ( force, dispDofs ) += elemForce1;
    select ( force, phiDofs  ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void PhaseFractureExtModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   state,
    const Vector&   state0,
    const Vector&   state00 )

{
  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of displacement and phase-field dofs

  const int   dispCount  = nodeCount * rank_;  
  const int   phiCount   = nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement and phase-field
  
  Vector      disp       ( dispCount );
  Vector      phi        ( phiCount );

  // old step element vector state:
  // phase-field
  
  Vector      phi0       ( phiCount );

  // old old step element vector state:
  // phase-field
  
  Vector      phi00       ( phiCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );
  Vector      stressP    ( strCount );

  // internal force vector:
  // displacement and phase-field

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( phiCount  );

  // element stiffness matrices (four components)

  Matrix      elemMat1   ( dispCount, dispCount ); // disp-disp
  Matrix      elemMat2   ( dispCount, phiCount  ); // disp-phi
  Matrix      elemMat3   ( phiCount,  dispCount ); // phi-disp
  Matrix      elemMat4   ( phiCount,  phiCount  ); // phi-phi
  
  // B matrices
  // displacement and phase-field
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, phiCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   phiDofs    ( phiCount  );
  IdxVector   dispDofs   ( dispCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  double      gphiEx;
  double      dgphi;
  double      ddgphi;

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
    dofs_->getDofIndices ( phiDofs, inodes, phiTypes_ );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the shape function

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current element displacements and phase-fields

    phi   = select ( state, phiDofs  );
    disp  = select ( state, dispDofs );

    // Get old step and old old step element phase-field

    phi0  = select ( state0,  phiDofs );
    phi00 = select ( state00, phiDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
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

      // Compute the displacement B-matrix for this integration point.

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // Get B-matrix associated with phase-field dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the strain for this integration point.

      matmul ( strain,  bd, disp  );

      // Store the regular strain components.

      strain_(slice(BEGIN,strCount),ipoint) = strain;

      // Compute the integration point phase-fields (current, old, old old)

      Vector Nip           = N( ALL,ip );

      double pf            = dot ( Nip, phi   );
      double pf0           = dot ( Nip, phi0  );
      double pf00          = dot ( Nip, phi00 );
      double pfEx          = 0.0;

      // Linearly extrapolate old old phase-field to current step

      if (dt0_ < 1e-20 )
      {
        pfEx = pf0;
      }
      else
      {
        // pfEx              = pf00 + ( dt_ + dt0_ )/dt0_ * ( pf0 - pf00 );
        // corrected Wick book
        pfEx = - pf00 * dt_/dt0_ + pf0 * (dt_+dt0_)/dt0_;
      }
      
      // Ensure extrapolated phase-field is within bounds [0,1]

      if ( pfEx < 0.0 )
      {
        pfEx = 0.0;
      }
      else if ( pfEx > 1.0 )
      {
        pfEx = 1.0;
      }

      // Compute the degradation function with the extrapolated phase-field

      {
        double gphiEx_n    = ::pow( 1- pfEx, p_);
        double gphiEx_d    = gphiEx_n + a1_*pfEx + a1_*a2_*pfEx*pfEx + 
                             a1_*a2_*a3_*pfEx*pfEx*pfEx;
        gphiEx             = gphiEx_n/gphiEx_d + 1.e-7;

      }

      strain_(strCount,ipoint) = gphiEx;

      // Compute the degradation function derivatives with the current phase-field

      {
        double gphi_n    = ::pow( 1- pf, p_);
        double gphi_d    = gphi_n + a1_*pf + a1_*a2_*pf*pf + a1_*a2_*a3_*pf*pf*pf;

        double dgphi_n   = - p_ * ( ::pow( 1- pf, p_- 1) );
        double dgphi_d   = dgphi_n + a1_ + 2.0*a1_*a2_*pf + 3.0*a1_*a2_*a3_*pf*pf;

        double ddgphi_n  = p_ * (p_ - 1) * ( ::pow( 1- pf, p_- 2) );
        double ddgphi_d  = ddgphi_n +  2.0*a1_*a2_ + 6.0*a1_*a2_*a3_*pf;

        dgphi            = ( dgphi_n*gphi_d - gphi_n*dgphi_d ) / ( gphi_d * gphi_d );
        ddgphi           = ( ( ddgphi_n * gphi_d - gphi_n * ddgphi_d ) * gphi_d - 2.0 * 
                           ( dgphi_n * gphi_d - gphi_n * dgphi_d ) * dgphi_d ) 
                         / ( gphi_d * gphi_d * gphi_d );  
      }

      // Compute the locally dissipated fracture energy function derivatives

      double dw        = eta_ + 2.0 * ( 1- eta_ ) * pf;
      double ddw       = 2.0 * ( 1- eta_ );

      // Compute stress and material tangent stiffness

      phaseMaterial_->update ( stress, stiff, strain, gphiEx, ipoint );

      // Get fracture driving energy and its derivative

      double Psi = phaseMaterial_->givePsi();
      stressP    = phaseMaterial_->giveDerivPsi();
      Psi        = max(Psi,Psi0_);

      // Compute stiffness matrix components

      wip         = ipWeights[ip];
      elemMat1   += wip * mc3.matmul ( bdt, stiff, bd );
      elemMat3   += wip * dgphi * ( matmul ( Nip, matmul(stressP, bd) ) );
      elemMat4   += wip * ( ( gc_/(cw_*l0_) * ddw + ddgphi * Psi) 
                    * matmul ( Nip, Nip ) + (2.0 * gc_*l0_/cw_) 
                    * mc2. matmul ( bet, be ) );
     
      // Compute internal forces

      elemForce1 +=  wip * ( mc1.matmul ( bdt, stress ) );
      elemForce2 +=  wip * ( Nip * ( gc_/(cw_*l0_) * dw + dgphi * Psi) 
                     + (2.0 * gc_*l0_/cw_) 
                     * mc1.matmul ( mc2. matmul ( bet, be ), phi ) );

      // Penalization (enforce fracture irreversibility) 

      if ( pf0 - pf > 1.e-15 )
      {
        double Dpf = pf0 - pf;
        elemMat4   += wip * ( gc_/(cw_*l0_) * beta_ * matmul ( Nip, Nip ) );
        elemForce2 -= wip * ( gc_/(cw_*l0_) * beta_ * Dpf ) * Nip;
      }   

    }  // End of loop on integration points

    // Assembly ...

    mbuilder.addBlock ( dispDofs, dispDofs, elemMat1 );
    mbuilder.addBlock ( dispDofs, phiDofs,  elemMat2 );
    mbuilder.addBlock ( phiDofs,  dispDofs, elemMat3 );
    mbuilder.addBlock ( phiDofs,  phiDofs,  elemMat4 );

    select ( force, dispDofs ) += elemForce1;
    select ( force, phiDofs  ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void PhaseFractureExtModel::getMatrix2_

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
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  PhaseFractureExtModel::initializeIPMPMap_ ( )

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


bool PhaseFractureExtModel::getTable_

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


void PhaseFractureExtModel::getStress_

  ( XTable&        table,
    const Vector&  weights )

{
  Matrix     sfuncs     = shape_->getShapeFunctions ();

  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];

  Matrix     stiff      ( strCount,  strCount );
  Matrix     ndStress   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  Vector     stressIp   ( strCount );
  Vector     strainIp   ( strCount );

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
      ipoint     = ipMpMap_ (ielem,ip);

      strainIp      = strain_(slice(BEGIN,strCount),igpoint);
      double gphiEx = strain_(strCount,igpoint);

      phaseMaterial_->update ( stressIp, stiff, strainIp, gphiEx, ipoint );

      // Extrapolate the integration point stresses to the nodes using
      // the transposed shape functions.

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


void PhaseFractureExtModel::getStrain_

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

      ndStrain  += matmul ( sfuncs(ALL,ip), strain_(slice(BEGIN,strCount),ipoint) );
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
//   getOutputData_
//-----------------------------------------------------------------------

/** all data at Gauss points are computed, however, only those requested
 *  in the vtk block in the .pro file will be written to vtu files.
 *  Example: vtk.data = "stress_yy | stress_xx" for xx and yy stress 
 *  components.
 */

void PhaseFractureExtModel::getOutputData_

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
      
      double gphi = strain_(strCount,igpoint);

      phaseMaterial_->update ( stress, stiff, strain, gphi, ipoint );
 
      //funcs = abs ( funcs );

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

void PhaseFractureExtModel::checkCommit_

  ( const Properties&  params )

{
  // System::info() << myName_ << " : check commit ... do nothing!\n";
}


//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void PhaseFractureExtModel::setStepSize_

  ( const Properties&  params )

{

  double       dt;
  double       dt0;

  params.get ( dt,  XProps::STEP_SIZE   );
  params.get ( dt0, XProps::STEP_SIZE_0 );

  dt_  = dt;
  dt0_ = dt0;
  
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newPhaseFractureExtModel
//-----------------------------------------------------------------------


static Ref<Model>     newPhaseFractureExtModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<PhaseFractureExtModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declarePhaseFractureExtModel
//-----------------------------------------------------------------------


void declarePhaseFractureExtModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "PhaseFractureExt", & newPhaseFractureExtModel );
}
