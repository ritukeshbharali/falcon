
/** @file SaturatedPorousFractureModel.cpp
 *  @brief Saturated porous media model with phase-field fracture. 
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
 *     - [24 May 2023] Added randomly distributed porosity,
 *       arc-length mode, phase-field specific BFGS mode 
 *       features (RB)  
 *     - [25 December 2023] removed getIntForce_,
 *       getMatrix_ returns the internal force if
 *       mbuilder = nullptr. Eliminates duplicate code. (RB)
 */

/* Include c++ headers */

#include <cstdlib> 

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/implict/ArclenActions.h>
#include <jive/geom/Geometries.h>
#include <jem/numeric/algebra/EigenUtils.h>
#include <jive/implict/SolverInfo.h>

/* Include other headers */

#include "FalconPoroMechModels.h"
#include "SaturatedPorousFractureModel.h"
#include "util/TensorUtils.h"
#include "util/TbFiller.h"
#include "util/XNames.h"

using jive::implict::ArclenActions;
using jive::implict::ArclenParams;

//=================================================================================================
//   class SaturatedPorousFractureModel
//=================================================================================================

//-------------------------------------------------------------------------------------------------
//   static data
//-------------------------------------------------------------------------------------------------


const char*  SaturatedPorousFractureModel::DISP_NAMES[3]          = { "dx", "dy", "dz" };
const char*  SaturatedPorousFractureModel::SHAPE_PROP             = "shape";
const char*  SaturatedPorousFractureModel::MATERIAL_PROP          = "material";

const char*  SaturatedPorousFractureModel::GENERATE_NEW_SEED      = "generateNewSeed";
const char*  SaturatedPorousFractureModel::KEEP_OFF_DIAGS_PROP    = "keepOffDiags";
const char*  SaturatedPorousFractureModel::ARCLEN_MODE_PROP       = "arcLenMode";
const char*  SaturatedPorousFractureModel::CONVEXIFY_PROP         = "convexify";

const char*  SaturatedPorousFractureModel::INTRIN_PERM_PROP       = "intrinPerm";
const char*  SaturatedPorousFractureModel::FLUID_VISC_PROP        = "fluidVisc";
const char*  SaturatedPorousFractureModel::SOLID_STIFF_PROP       = "solidStiff";
const char*  SaturatedPorousFractureModel::FLUID_STIFF_PROP       = "fluidStiff";
const char*  SaturatedPorousFractureModel::POROSITY_PROP          = "porosity";
const char*  SaturatedPorousFractureModel::BIOT_COEFF_PROP        = "biotCoeff";
const char*  SaturatedPorousFractureModel::DTIME_PROP             = "dtime";

const char*  SaturatedPorousFractureModel::SOLID_DENSITY_PROP     = "rhoSolid";
const char*  SaturatedPorousFractureModel::FLUID_DENSITY_PROP     = "rhoFluid";
const char*  SaturatedPorousFractureModel::GRAVITY_PROP           = "gravity";

const char*  SaturatedPorousFractureModel::FRACTURE_POROSITY_PROP = "fracPorosity";
const char*  SaturatedPorousFractureModel::RANDOM_POROSITY_PROP   = "randPorosity";

const char*  SaturatedPorousFractureModel::PERM_TYPE_PROP         = "permType";
const char*  SaturatedPorousFractureModel::FIXED_ISO_PERM         = "fixedIsoPerm";
const char*  SaturatedPorousFractureModel::VOL_STRAIN_PERM        = "volStrainPerm";
const char*  SaturatedPorousFractureModel::TENSILE_STRAIN_PERM    = "tensileStrainPerm";
const char*  SaturatedPorousFractureModel::CUBIC_PERM             = "cubicPerm";
const char*  SaturatedPorousFractureModel::CUBIC_ISO_PERM         = "cubicIsoPerm";

const char*  SaturatedPorousFractureModel::PRESSURE_PSI_PROP      = "pressurePsi";
const char*  SaturatedPorousFractureModel::FRACTURE_TYPE_PROP     = "fractureType";
const char*  SaturatedPorousFractureModel::GRIFFITH_ENERGY_PROP   = "griffithEnergy";
const char*  SaturatedPorousFractureModel::LENGTH_SCALE_PROP      = "lengthScale";
const char*  SaturatedPorousFractureModel::TENSILE_STRENGTH_PROP  = "tensileStrength";

const char*  SaturatedPorousFractureModel::BRITTLE_AT1             = "at1";
const char*  SaturatedPorousFractureModel::BRITTLE_AT2             = "at2";
const char*  SaturatedPorousFractureModel::LINEAR_CZM              = "czmLinear";
const char*  SaturatedPorousFractureModel::EXPONENTIAL_CZM         = "czmExponential";
const char*  SaturatedPorousFractureModel::CORNELISSEN_CZM         = "czmCornelissen";
const char*  SaturatedPorousFractureModel::HYPERBOLIC_CZM          = "czmHyperbolic";

vector<String> SaturatedPorousFractureModel::fractureModels      (6);
vector<String> SaturatedPorousFractureModel::permeabilityModels  (5);

//-------------------------------------------------------------------------------------------------
//   constructors & destructor
//-------------------------------------------------------------------------------------------------


SaturatedPorousFractureModel::SaturatedPorousFractureModel

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


  // Get the DOF space, add displacement, pressure, phase-field DOFs

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( rank_ + 2 );

  // Make a shallow copy of dofTypes_ for displacement.

  dispTypes_.ref ( dofTypes_[slice(BEGIN,rank_)] );

  for ( int i = 0; i < rank_; i++ )
  {
    dispTypes_[i] = dofs_->addType ( DISP_NAMES[i] );
  }

  // Make a shallow copy of dofTypes_ for pressure.

  presTypes_.ref ( dofTypes_[slice(rank_,rank_+1)] );

  presTypes_[0] = dofs_->addType ( "dp" );

  // Make a shallow copy of dofTypes_ for phase-field.

  phiTypes_.ref ( dofTypes_[slice(rank_+1,END)] );

  phiTypes_[0] = dofs_->addType ( "phi" );
    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  // Allocate memory for the strains and set them to zero

  strain_.resize ( STRAIN_COUNTS[rank_] + 2, ipCount );
  strain_ = 0.0;

  // Create a material model object.

  material_      = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );
  phaseMaterial_ = dynamicCast<PhaseFractureMaterial> ( material_ );
  material_-> allocPoints  ( ipCount );

  initializeIPMPMap_ ();

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

  n0_         = 0.375;
  myProps.find ( n0_, POROSITY_PROP );
  myConf. set  ( POROSITY_PROP, n0_ );

  alpha_      = 1.0;
  myProps.find ( alpha_, BIOT_COEFF_PROP );
  myConf. set  ( BIOT_COEFF_PROP, alpha_ );

  dtime_       = 1.0;
  myProps.find ( dtime_, DTIME_PROP );
  myConf. set  ( DTIME_PROP, dtime_ );

  fracPorosity_ = false;
  myProps.find ( fracPorosity_, FRACTURE_POROSITY_PROP );
  myConf. set  ( FRACTURE_POROSITY_PROP, fracPorosity_ );

  rhoS_        = 2000.0;
  myProps.find ( rhoS_, SOLID_DENSITY_PROP );
  myConf. set  ( SOLID_DENSITY_PROP, rhoS_ );

  rhoF_        = 1000.0;
  myProps.find ( rhoF_, FLUID_DENSITY_PROP );
  myConf. set  ( FLUID_DENSITY_PROP, rhoF_ );

  //  Setup the permeability models

  permeabilityModels[0]     = FIXED_ISO_PERM;
  permeabilityModels[1]     = VOL_STRAIN_PERM;
  permeabilityModels[2]     = TENSILE_STRAIN_PERM;
  permeabilityModels[3]     = CUBIC_PERM;
  permeabilityModels[4]     = CUBIC_ISO_PERM;

  // Get permeability type

  myProps.get ( permType_, PERM_TYPE_PROP  );

  if ( std::find ( permeabilityModels.begin (),
                   permeabilityModels.end   (),
                   permType_ ) == permeabilityModels.end () )
  {
    throw Error (
      JEM_FUNC,
      String("unexpected definition of permeability type!!!\n") +
      String("Supported types include: \n") +
      FIXED_ISO_PERM   + String(", ") + VOL_STRAIN_PERM + String(", ") +
      TENSILE_STRAIN_PERM + String(", ") + CUBIC_PERM + String(", ") +
      CUBIC_ISO_PERM
    );
  }

  myConf. set  ( PERM_TYPE_PROP, permType_ );

  bool generateNewSeed = false;
  myProps.find ( generateNewSeed, GENERATE_NEW_SEED );

  if ( generateNewSeed )
  {
    // Initialize the seed with time(0), the current unix timestamp
    // On every run, a new set of random numbers are generated.

    srand(time(0));
  }

  // Fill in random values (0.0,1.0]

  randPorosity_ = true;
  nrand_.resize( ipCount );

  myProps.find ( randPorosity_, RANDOM_POROSITY_PROP );
  myConf. set  ( RANDOM_POROSITY_PROP, randPorosity_ );

  if ( randPorosity_ )
  {
    std::generate(nrand_.begin(),nrand_.end(), rand );
    nrand_ /= max(nrand_);
    nrand_ -= 0.5;
    nrand_ *= 2.0;
  }
  else
  {
    nrand_ = 0.0;
  }

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

  // Keep diagonal terms of the stiffness matrix

  keepOffDiags_  = false;
  myProps.find ( keepOffDiags_, KEEP_OFF_DIAGS_PROP );
  myConf. set  ( KEEP_OFF_DIAGS_PROP, keepOffDiags_ );

  // Model contribution to arc-length function

  arcLenMode_  = false;
  myProps.find ( arcLenMode_, ARCLEN_MODE_PROP );
  myConf. set  ( ARCLEN_MODE_PROP, arcLenMode_ );

  // Setup the fracture models

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

  pressurePsi_ = false;
  myProps.find ( pressurePsi_, PRESSURE_PSI_PROP );
  myConf. set  ( PRESSURE_PSI_PROP, pressurePsi_ );

  gc_          = 2.7;
  myProps.find ( gc_, GRIFFITH_ENERGY_PROP );
  myConf. set  ( GRIFFITH_ENERGY_PROP, gc_ );

  l0_          = 0.015;
  myProps.find ( l0_, LENGTH_SCALE_PROP );
  myConf. set  ( LENGTH_SCALE_PROP, l0_ );

  sigmaT_      = 0.0;
  myProps.find ( sigmaT_, TENSILE_STRENGTH_PROP );
  myConf. set  ( TENSILE_STRENGTH_PROP, sigmaT_ );

  convexify_ = false;
  myProps.find ( convexify_, CONVEXIFY_PROP );
  myConf. set  ( CONVEXIFY_PROP, convexify_ );
  
  preHist_.H0_.resize ( ipCount );
  newHist_.H0_.resize ( ipCount );
 
}


SaturatedPorousFractureModel::~SaturatedPorousFractureModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void SaturatedPorousFractureModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  material_-> configure ( matProps, globdat );

  // Compute derived quantities

  const double E    = material_->getYoung   ();

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
  else if ( fracType_ == LINEAR_CZM )
  {
    cw_    = 3.14159265359;
    eta_   = 2.0;
    p_     = 2.0;
    a1_    = ( 4.0 * E * gc_ )/( cw_ * l0_ * sigmaT_ * sigmaT_ );
    a2_    = -0.5;
    a3_    = 0.0;
    Psi0_  = 0.5 * (sigmaT_ * sigmaT_)/E;
  }
  else if ( fracType_ == EXPONENTIAL_CZM )
  {
    cw_    = 3.14159265359;
    eta_   = 2.0;
    p_     = 2.5;
    a1_    = ( 4.0 * E * gc_ )/( cw_ * l0_ * sigmaT_ * sigmaT_ );
    a2_    = ::pow(2.0, 5.0/3.0) - 3.0;
    a3_    = 0.0;
    Psi0_  = 0.5 * (sigmaT_ * sigmaT_)/E;
  }
  else if ( fracType_ == CORNELISSEN_CZM )
  {
    cw_    = 3.14159265359;
    eta_   = 2.0;
    p_     = 2.0;
    a1_    = ( 4.0 * E * gc_ )/( cw_ * l0_ * sigmaT_ * sigmaT_ );
    a2_    = 1.3868;
    a3_    = 0.6567;
    Psi0_  = 0.5 * (sigmaT_ * sigmaT_)/E;
  }
  else if ( fracType_ == HYPERBOLIC_CZM )
  {
    cw_    = 3.14159265359;
    eta_   = 2.0;
    p_     = 4.0;
    a1_    = ( 4.0 * E * gc_ )/( cw_ * l0_ * sigmaT_ * sigmaT_ );
    a2_    = ::pow(2.0, 7.0/3.0) - 4.5;
    a3_    = 0.0;
    Psi0_  = 0.5 * (sigmaT_ * sigmaT_)/E;
  }

  // Setup history variables

  preHist_.H0_  = Psi0_;
  newHist_.H0_  = 0.0;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void SaturatedPorousFractureModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP  );

  material_-> getConfig ( matConf, globdat );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool SaturatedPorousFractureModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_INT_VECTOR )
  {
    Vector  state;
    Vector  state0;
    Vector  force;

    // Get the current and old step displacements.

    StateVector::get    ( state, dofs_, globdat  );
    StateVector::getOld ( state0, dofs_, globdat );

    // Get the internal force vector.

    params.get ( force, ActionParams::INT_VECTOR );

    getMatrix_ ( nullptr, force, state, state0 );

    globdat.set ( XProps::FE_DISSIPATION, diss_ );

    return true;
  }

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {

    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  state0;
    Vector  force;

    // Get the current and old step displacements.

    StateVector::get    ( state, dofs_, globdat  );
    StateVector::getOld ( state0, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( mbuilder, force, state, state0 );

    globdat.set ( XProps::FE_DISSIPATION, diss_ );

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

  /*  During initialization of the model, the pressure dofs of the
   *  mid nodes of all domain elements are linearly constrained to
   *  the corner nodes of the respective edge. This yields a Taylor
   *  Hood element.
  */ 

  if ( action == Actions::INIT )
  {
    // Get the constraints associated with the DOF space.

    Ref<Constraints>  cons = Constraints::get ( dofs_, globdat );

    setTaylorHoodConstraints_ ( *cons );

    globdat.set ( Globdat::TIME, 0.0 );

    return true;
  }

  /*  COMMIT: Requests an action when a computation step has converged.
   *  At this point, the material internal state variables are updated.
   */ 

  if ( action == Actions::COMMIT )
  {
    // material_->commit ();

    newHist_.H0_. swap ( preHist_.H0_  );

    double t;

    globdat.get ( t, Globdat::TIME );

    t += dtime_;

    globdat.set ( Globdat::TIME, t );


    return true;
  }

  /*  CHECK_COMMIT: Can used to discard the current step
  */ 

  if ( action == Actions::CHECK_COMMIT )
  {
    // checkCommit_ ( params );
    return true;
  }

  /*  GET_TABLE: Requests a post-processing action, such as computing
   *  tabular data (stress, strains, etc.) to print VTK files.
   */ 
  
  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  /*  SET_STEP_SIZE: Implements changes to the step-size, thrown by a
   *  solver module (e.g., AdaptiveStepping, ReduceStepping)
   */ 

  if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  // GET_ARC_FUNC: Requests the arc-length function
  
  if ( action == ArclenActions::GET_ARC_FUNC && arcLenMode_ )
  {
    getArcFunc_ ( params, globdat );

    return true;
  }

  // GET_UNIT_LOAD: Requests the unit load
  
  if ( action == ArclenActions::GET_UNIT_LOAD && arcLenMode_ )
  {
    getUnitLoad_ ( params, globdat );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void SaturatedPorousFractureModel::getMatrix_

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
  // displacement, pressure, phase-field
  
  Vector      disp       ( dispCount );
  Vector      pres       ( nodeCount );
  Vector      phi        ( nodeCount );

  // old step element vector state:
  // displacement, pressure, phase-field
  
  Vector      disp0      ( dispCount );
  Vector      pres0      ( nodeCount );
  Vector      phi0       ( nodeCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      stressP    ( strCount );
  Vector      strain     ( strCount );
  Vector      strain0    ( strCount );
  

  // internal force vector:
  // displacement and pressure

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( nodeCount  );
  Vector      elemForce3 ( nodeCount );

  // element stiffness matrices (nine components)

  Matrix      elemMat11  ( dispCount, dispCount ); // disp-disp
  Matrix      elemMat12  ( dispCount, nodeCount ); // disp-pres
  Matrix      elemMat13  ( dispCount, nodeCount ); // disp-phi

  Matrix      elemMat21  ( nodeCount, dispCount ); // pres-disp
  Matrix      elemMat22  ( nodeCount, nodeCount ); // pres-pres
  Matrix      elemMat23  ( nodeCount, nodeCount ); // pres-phi

  Matrix      elemMat31  ( nodeCount, dispCount );  // phi-disp
  Matrix      elemMat32  ( nodeCount, nodeCount );  // phi-pres
  Matrix      elemMat33  ( nodeCount, nodeCount );  // phi-phi
  
  // B matrices
  // displacement and pressure
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, nodeCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   dispDofs   ( dispCount );
  IdxVector   presDofs   ( nodeCount );
  IdxVector   phiDofs    ( nodeCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Initialize degradation function and its derivatives

  double      gphi;
  double      gphi0;
  double      dgphi;
  double      ddgphi;

  // Initialize dissipation

  diss_ = 0.0;

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
    dofs_->getDofIndices ( phiDofs,  inodes, phiTypes_  );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();

    // Compute the area of the element

    const double area = sum( ipWeights );
    
    // Get current nodal displacements and pressures

    disp = select ( state, dispDofs );
    pres = select ( state, presDofs );
    phi  = select ( state, phiDofs  );

    // Get old step nodal displacements and pressures

    disp0 = select ( state0, dispDofs );
    pres0 = select ( state0, presDofs );
    phi0  = select ( state0,  phiDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0;
    elemForce3 = 0.0;

    if ( mbuilder != nullptr )
    {
      elemMat11  = 0.0;
      elemMat12  = 0.0;
      elemMat13  = 0.0;

      elemMat21  = 0.0;
      elemMat22  = 0.0;
      elemMat23  = 0.0;

      elemMat31  = 0.0;
      elemMat32  = 0.0;
      elemMat33  = 0.0;
    }

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
      matmul ( strain0, bd, disp0 );

      const double evol  = dot ( voigtDiv_, strain  );
      const double evol0 = dot ( voigtDiv_, strain0 );

      // Store the regular strain components.

      strain_(slice(BEGIN,strCount),ipoint) = strain;

      // Compute the integration point pressures

      const Vector Nip  = N( ALL,ip );
      const double p    = dot ( Nip, pres  );
      const double p0   = dot ( Nip, pres0 );

      // Compute the integration point phase-fields (current, old, old old)

      double pf  = dot ( Nip, phi  );
      double pf0 = dot ( Nip, phi0 );

      // Compute the degradation function with the extrapolated phase-field

      {
        double gphi0_n    = ::pow( 1.0 - pf0, p_);
        double gphi0_d    = gphi0_n + a1_*pf0 + a1_*a2_*pf0*pf0 + 
                             a1_*a2_*a3_*pf0*pf0*pf0;
        gphi0             = gphi0_n/gphi0_d + 1.e-6;
      }

      // Compute the degradation function derivatives with the current phase-field

      {
        double gphi_n    = ::pow( 1.0- pf, p_);
        double gphi_d    = gphi_n + a1_*pf + a1_*a2_*pf*pf + a1_*a2_*a3_*pf*pf*pf;

        double dgphi_n   = - p_ * ( ::pow( 1.0- pf, p_- 1.0) );
        double dgphi_d   = dgphi_n + a1_ + 2.0*a1_*a2_*pf + 3.0*a1_*a2_*a3_*pf*pf;

        double ddgphi_n  = p_ * (p_ - 1.0) * ( ::pow( 1.0- pf, p_- 2.0) );
        double ddgphi_d  = ddgphi_n +  2.0*a1_*a2_ + 6.0*a1_*a2_*a3_*pf;

        gphi             = gphi_n/gphi_d + 1.e-6;
        dgphi            = ( dgphi_n*gphi_d - gphi_n*dgphi_d ) / ( gphi_d * gphi_d );
        ddgphi           = ( ( ddgphi_n * gphi_d - gphi_n * ddgphi_d ) * gphi_d - 2.0 * 
                           ( dgphi_n * gphi_d - gphi_n * dgphi_d ) * dgphi_d ) / ( gphi_d * gphi_d * gphi_d );  
      }

      strain_(strCount, ipoint) = gphi;

      // Compute the locally dissipated fracture energy function derivatives

      double dw  = eta_ + 2.0 * ( 1- eta_ ) * pf;
      double ddw = 2.0 * ( 1- eta_ );

      // Compute stress and material tangent stiffness

      if ( convexify_ )
      {
        phaseMaterial_->update ( stress, stiff, strain, gphi0, ipoint );
      }
      else
      {
        phaseMaterial_->update ( stress, stiff, strain, gphi, ipoint );
      }

      // Get the fracture driving energy and its derivative

      double Psi = phaseMaterial_->givePsi();
      stressP    = phaseMaterial_->giveDerivPsi();

      // Compute the storage coefficient and its derivative

      const double Sto = ( ( alpha_ - n0_ ) / Ks_ ) + ( n0_ / Kf_ );

      // Compute dissipation

      diss_ += ipWeights[ip] * ( gphi0 - gphi ) * Psi;   // CHANGED

      // Get old history value and update fracture driving energy

      const double H0      = preHist_.H0_[ipoint];
      double       loading = 1.0;

      if ( Psi < H0 || Psi < Psi0_ )
      {
        loading = 0.0;
      }

      Psi                  = max(Psi,H0);
      Psi                  = max(Psi,Psi0_);
      newHist_.H0_[ipoint] = Psi;

      // Compute bulk and fracture permeability

      Matrix kappaComb ( rank_, rank_ );  kappaComb = 0.0;
      Matrix kappaBulk ( rank_, rank_ );  kappaBulk = 0.0;
      Matrix kappaFrac ( rank_, rank_ );  kappaFrac = 0.0;
      Matrix I2        ( rank_, rank_ );  I2        = 0.0;

      I2         = tensorUtils::getI2 ( rank_ );

      // Permeability models

      if ( permType_ == FIXED_ISO_PERM )
      {
        kappaComb  = Keff_ * I2;
      }
      else if ( permType_ == VOL_STRAIN_PERM )
      {
        double effStrain = 0.3333 * ( strain[0] + strain[1] + strain[2] ) ;

        if ( effStrain > 0.0 )
        {
          kappaComb = ( Keff_ + 1.e+8 * effStrain * Keff_ ) * I2;
        }
        else
        {
          kappaComb  = Keff_ * I2;
        }
      }
      else if ( permType_ == TENSILE_STRAIN_PERM )
      {
        Matrix strainT = tensorUtils::voigt2tensorStrain( strain );

        Matrix  eigVecs;    eigVecs .resize ( 3, 3 );
        Vector  eigVals;    eigVals .resize ( 3 );

        using jem::numeric::EigenUtils;

        EigenUtils::symSolve ( eigVals, eigVecs, strainT );

        for ( int is = 0; is < 3 ; is++ )
        {
          if ( eigVals[is] < 0.0 )
          {
            eigVals[is] = 0.0;
          }
        }

        double effStrain = jem::min ( eigVals );

        if ( effStrain > 0.0 )
        {
          kappaComb = ( Keff_ + 1.e+7 * effStrain * Keff_ ) * I2;
        }
        else
        {
          kappaComb  = Keff_ * I2;
        }

      }
      else // cubic (aniso and iso versions)
      {
        kappaBulk = Keff_ * I2;

        // Compute residual crack opening 

        double wr = ::sqrt( 120. * kappa_ ); // sqrt ( 12 * 10 * kappa_ )

        // Compute normalized phase-field gradient

        Vector   pfGrad  ( rank_ ); pfGrad = 0.0;
        matmul ( pfGrad, be, phi );

        double  pfGradNorm = jem::numeric::norm2(pfGrad);
        pfGrad /= pfGradNorm;

        // Get strain in a tensorial format

        Matrix strainT = tensorUtils::voigt2tensorRankStrain ( strain );

        // Compute the fracture opening

        double wc = ::sqrt(area) * ( 1.0 + dot( pfGrad, matmul( strainT, pfGrad ) ) );

        // Compute final fracture opening ( max of wc and wr )

        double wh = max( wr, wc );

        // Compute fracture permeability

        kappaFrac = ( wh * wh / ( 12. * mu_ ) ) * I2 ;

        // Add anisotropy

        if ( permType_ == CUBIC_PERM )
        {
          kappaFrac -= ( wh * wh / ( 12. * mu_ ) ) * matmul ( pfGrad, pfGrad ) ;
        }

        // Compute dual permeability

        if ( pf > 0.8 )
        {
          kappaComb = kappaBulk + 25. * ::pow( pf - 0.8, 2.0 ) * kappaFrac;
        }
        else
        {
          kappaComb = kappaBulk;
        }
      }

      // Compute weight of the ip

      wip         = ipWeights[ip];

      // Compute stiffness matrix components

      if ( mbuilder != nullptr )
      {
        elemMat11  += wip * mc3.matmul ( bdt, stiff, bd );
        elemMat12  -= wip * alpha_ * mc2.matmul ( bdt, matmul (voigtDiv_, Nip ) );

        if ( !convexify_ )
        {
          elemMat13  += wip * dgphi * ( matmul ( matmul(bdt, stressP), Nip ) );
        }

        elemMat21  += wip * alpha_ / dtime_ * mc2.matmul ( matmul (Nip, voigtDiv_), bd );
        elemMat22  += wip * (  Sto / dtime_ * matmul ( Nip, Nip )
                             + mc3.matmul ( bet, kappaComb, be ) );

        elemMat31  += wip * loading * dgphi * ( matmul ( Nip, matmul(stressP, bd) ) );
        elemMat33  += wip * ( ( gc_/(cw_*l0_) * ddw + ddgphi * Psi) 
                          * matmul ( Nip, Nip ) + (2.0 * gc_*l0_/cw_) 
                          * mc2. matmul ( bet, be ) );
      }

      // compute internal forces

      elemForce1 +=  wip * ( mc1.matmul ( bdt, stress ) );
      elemForce1 -=  wip * alpha_ * p * ( mc1.matmul ( bdt, voigtDiv_ ) );
      elemForce2 +=  wip *  (  Nip * (  Sto * (p-p0) + alpha_ * ( evol - evol0 )  ) / dtime_ 
                             + mc1.matmul ( mc3.matmul( bet, kappaComb, be ), pres ) );
      elemForce3 +=  wip * ( Nip * ( gc_/(cw_*l0_) * dw + dgphi * Psi) + (2.0 * gc_*l0_/cw_) * mc1.matmul ( mc2. matmul ( bet, be ), phi ) );

    }  // end of loop on integration points

    // Assembly ...

    if ( mbuilder != nullptr )
    {
      mbuilder -> addBlock ( dispDofs, dispDofs, elemMat11 );
      mbuilder -> addBlock ( dispDofs, presDofs, elemMat12 );
      mbuilder -> addBlock ( presDofs, dispDofs, elemMat21 );
      mbuilder -> addBlock ( presDofs, presDofs, elemMat22 );
      mbuilder -> addBlock ( phiDofs,  phiDofs,  elemMat33 );

      // Assemble the off-diagonal components involving the phase-field
      // For staggered solvers, keepOffDiags_ = false;

      if ( keepOffDiags_ )
      {
        mbuilder -> addBlock ( dispDofs, phiDofs,  elemMat13 );
        mbuilder -> addBlock ( phiDofs,  dispDofs, elemMat31 );
        mbuilder -> addBlock ( presDofs, phiDofs,  elemMat23 );
        mbuilder -> addBlock ( phiDofs,  presDofs, elemMat32 );  
      }
    }

    select ( force, dispDofs ) += elemForce1;
    select ( force, presDofs ) += elemForce2;
    select ( force, phiDofs  ) += elemForce3;
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void SaturatedPorousFractureModel::getMatrix2_

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
      // compute matrix of shape function N       

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Assemble mass matrix

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}

//-----------------------------------------------------------------------
//   getArcFunc
//-----------------------------------------------------------------------

void SaturatedPorousFractureModel::getArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::StateVector;
  using jive::implict::ArclenActions;
  using jive::implict::ArclenParams;

  // Get the state vectors (current and old)

  Vector  state;
  Vector  state0;

  StateVector::get       ( state,   dofs_, globdat );
  StateVector::getOld    ( state0,  dofs_, globdat );

  // Get the arc-length parameters and set them to zero

  Vector  jac10; 
  double  jac11;

  params.get ( jac10,   ArclenParams::JACOBIAN10 );

  jac10 = 0.0;
  jac11 = 0.0;

  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of displacement and phase-field dofs

  const int   dispCount  = nodeCount * rank_;  

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement and phase-field
  
  Vector      disp       ( dispCount );
  Vector      phi        ( nodeCount );

  Vector      phi0       ( nodeCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );
  Vector      stressP    ( strCount );

  // internal force vector:
  // displacement, pressure, phase-field

  Vector      elemjac10_0 ( dispCount );
  Vector      elemjac10_1 ( nodeCount );
  Vector      elemjac10_2 ( nodeCount );
  
  // B matrices
  // displacement and phase-field
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, nodeCount);
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   phiDofs    ( nodeCount );
  IdxVector   dispDofs   ( dispCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  double      gphi;
  double      gphi0;
  double      dgphi;

  double      fvalue = 0.0;

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
    dofs_->getDofIndices ( phiDofs,  inodes, phiTypes_  );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Get the shape function

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current nodal displacements and pressures

    disp = select ( state, dispDofs );
    phi  = select ( state, phiDofs  );

    // Get old step nodal displacements and pressures

    phi0  = select ( state0,  phiDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemjac10_0 = 0.0;
    elemjac10_1 = 0.0;
    elemjac10_2 = 0.0;   

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

      // Compute the integration point phase-fields (current, old, old old)

      const Vector Nip  = N( ALL,ip );
      double pf         = dot ( Nip, phi  );
      double pf0        = dot ( Nip, phi0 );

      // Compute the degradation function with the extrapolated phase-field

      {
        pf0 = max ( pf0, 0.0 );
        pf0 = min ( pf0, 1.0 );

        double gphi0_n    = ::pow( 1.0 - pf0, p_);
        double gphi0_d    = gphi0_n + a1_*pf0 + a1_*a2_*pf0*pf0 + 
                             a1_*a2_*a3_*pf0*pf0*pf0;
        gphi0             = gphi0_n/gphi0_d + 1.e-6;
      }

      // Compute the degradation function derivatives with the current phase-field

      {
        pf = max ( pf, 0.0 );
        pf = min ( pf, 1.0 );

        double gphi_n    = ::pow( 1.0- pf, p_);
        double gphi_d    = gphi_n + a1_*pf + a1_*a2_*pf*pf + a1_*a2_*a3_*pf*pf*pf;

        double dgphi_n   = - p_ * ( ::pow( 1.0- pf, p_- 1.0) );
        double dgphi_d   = dgphi_n + a1_ + 2.0*a1_*a2_*pf + 3.0*a1_*a2_*a3_*pf*pf;

        gphi             = gphi_n/gphi_d + 1.e-6;
        dgphi            = ( dgphi_n*gphi_d - gphi_n*dgphi_d ) / ( gphi_d * gphi_d );
      }

      // Compute stress and material tangent stiffness

      phaseMaterial_->update ( stress, stiff, strain, gphi, ipoint );

      // Get the fracture driving energy and its derivative

      double Psi = phaseMaterial_->givePsi();
      stressP    = phaseMaterial_->giveDerivPsi();

      // Compute stiffness matrix components

      wip         = ipWeights[ip];

      fvalue      +=  wip * ( gphi0 - gphi ) * Psi;
      elemjac10_0 +=  wip * ( gphi0 - gphi ) * ( mc1.matmul ( bdt, stressP ) );
      elemjac10_2 -=  wip * Nip * dgphi * Psi;

    }  // End of loop on integration points

    // Assembly ...

    select ( jac10, dispDofs ) += elemjac10_0;
    select ( jac10, phiDofs  ) += elemjac10_2;

  }

  globdat.set ( XProps::FE_DISSIPATION,    fvalue );
  params .set ( ArclenParams::JACOBIAN10,  jac10  );
  params .set ( ArclenParams::JACOBIAN11,  jac11  );

}

//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------

void SaturatedPorousFractureModel::getUnitLoad_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::StateVector;
  using jive::implict::ArclenActions;
  using jive::implict::ArclenParams;

  idx_t       iiter = 0; 
  globdat.get (  iiter , "iiter" );

  // Get the state vectors (current and old)

  Vector  state;
  Vector  state0;

  StateVector::get       ( state,   dofs_, globdat );

  if ( iiter == 0 )
  {
    StateVector::getOldOld ( state0,  dofs_, globdat );
  }
  else
  {
    StateVector::getOld    ( state0,  dofs_, globdat );
  }
  

  // Get the arc-length parameters and set them to zero

  Vector  unitLoad; 

  params.get ( unitLoad, ArclenParams::UNIT_LOAD );

  unitLoad = 0.0;

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
  // displacement, pressure, phase-field
  
  Vector      disp       ( dispCount );
  Vector      pres       ( nodeCount );
  //Vector      phi        ( nodeCount );

  // old step element vector state:
  // displacement, pressure, phase-field
  
  Vector      disp0      ( dispCount );
  Vector      pres0      ( nodeCount );
  //Vector      phi0       ( nodeCount );

  // current and old strains

  Vector      strain     ( strCount );
  Vector      strain0    ( strCount );

  // internal force vector:
  // displacement and pressure

  Vector      elemUnitForce1 ( dispCount );
  Vector      elemUnitForce2 ( nodeCount );
  Vector      elemUnitForce3 ( nodeCount );
  
  // B matrices
  // displacement and pressure
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, nodeCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   dispDofs   ( dispCount );
  IdxVector   presDofs   ( nodeCount );
  IdxVector   phiDofs    ( nodeCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Initialize degradation function

  //double      gphi;

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
    dofs_->getDofIndices ( phiDofs,  inodes, phiTypes_  );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current nodal displacements and pressures

    disp = select ( state, dispDofs );
    pres = select ( state, presDofs );
    //phi  = select ( state, phiDofs  );

    // Get old step nodal displacements and pressures

    disp0 = select ( state0, dispDofs );
    pres0 = select ( state0, presDofs );
    //phi0  = select ( state0,  phiDofs  );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemUnitForce1 = 0.0;
    elemUnitForce2 = 0.0;
    elemUnitForce3 = 0.0; 

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
      matmul ( strain0, bd, disp0 );

      // Store the regular strain components.

      strain_(slice(BEGIN,strCount),ipoint) = strain;

      // Compute the integration point pressures

      const Vector Nip  = N( ALL,ip );
      const double p    = dot ( Nip, pres  );
      const double p0   = dot ( Nip, pres0 );

      // Compute the storage coefficients

      const double Sto   = ( ( alpha_ - n0_ ) / Ks_ ) + ( n0_ / Kf_ );

      // Compute volumetric strain

      double evol  = dot( voigtDiv_, strain  );
      double evol0 = dot( voigtDiv_, strain0 );

      // compute stiffness matrix components   

      wip             = ipWeights[ip];

      elemUnitForce2 -= wip * ( Nip * (1.0/(dtime_*dtime_)) * 
                              ( Sto * ( p - p0 ) + alpha_ * ( evol - evol0 ) ) );

    }  // end of loop on integration points

    // Assembly ...

    select ( unitLoad, dispDofs ) += elemUnitForce1;
    select ( unitLoad, presDofs ) += elemUnitForce2;
    select ( unitLoad, phiDofs  ) += elemUnitForce3;
  }

  // double norm2UnitLoad = jem::numeric::norm2( unitLoad );

  //System::out() << "Norm of Unit Load = " << norm2UnitLoad << "\n";

  globdat.set ( XProps::FE_UNIT_LOAD, unitLoad );

}

//-----------------------------------------------------------------------
//   setTaylorHoodConstraints_
//-----------------------------------------------------------------------


void SaturatedPorousFractureModel::setTaylorHoodConstraints_ 

    ( Constraints& cons )
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


void  SaturatedPorousFractureModel::initializeIPMPMap_ ( )

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


bool SaturatedPorousFractureModel::getTable_

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


void SaturatedPorousFractureModel::getStress_

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


void SaturatedPorousFractureModel::getStrain_

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


void SaturatedPorousFractureModel::getHistory_

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

void SaturatedPorousFractureModel::getOutputData_

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

void SaturatedPorousFractureModel::checkCommit_

  ( const Properties&  params )

{
  System::info() << myName_ << " : check commit ... do nothing!\n";
 
}

//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void SaturatedPorousFractureModel::setStepSize_

  ( const Properties&  params )

{
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
//   newSaturatedPorousFractureModel
//-----------------------------------------------------------------------


static Ref<Model>     newSaturatedPorousFractureModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<SaturatedPorousFractureModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSaturatedPorousFractureModel
//-----------------------------------------------------------------------


void declareSaturatedPorousFractureModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "SaturatedPorousFracture", & newSaturatedPorousFractureModel );
}
