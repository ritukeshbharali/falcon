
/** @file SaturatedPorousMicroFractureModel.h
 *  @brief Implements saturated porous media model with fracture.
 *  
 *  This class implements a mass conserving saturated two-phase
 *  porous media model together with phase-field fracture. Model
 *  options include randomly distributed porosity, different
 *  permeability models, arc-length solver mode ( assembles unit
 *  load, arc-length constraint ), BFGS mode ( skip phase-field 
 *  coupled terms in the stiffness matrix, keepOffDiags = false ).
 *  The model also updates the time in globdat, so that Paraview
 *  Module can print time stamps to a pvd file.
 * 
 *  Please consider citing the corresponding article
 *  (doi:10.1007/s00466-023-02380-1), if the code benefits you.
 *   
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
 */

/* Include c++ headers */

#include <vector>

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/base/Error.h>
#include <jem/base/System.h>
#include <jem/base/IllegalInputException.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/sparse/utilities.h>
#include <jem/util/Properties.h>
#include <jem/util/Flex.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/geom/Geometries.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/StdSquare.h>
#include <jive/geom/StdCube.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/geom/StdShapeFactory.h>
#include <jive/geom/StdBezierShape.h>
#include <jive/geom/ParametricArea.h>
#include <jive/geom/ParametricVolume.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Globdat.h>

/* Include falcon headers */

#include "util/BasicUtils.h"
#include "materials/Material.h"
#include "materials/PhaseFractureMaterial.h"

using namespace jem;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::SparseArray;
using jem::util::Properties;
using jem::util::Flex;
using jive::Vector;
using jive::IdxVector;
using jive::IntMatrix;
using jive::Matrix;
using jive::Cubix;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::util::Globdat;
using std::vector;


//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef MatmulChain<double,1>  MChain1;
typedef MatmulChain<double,2>  MChain2;
typedef MatmulChain<double,3>  MChain3;
typedef MatmulChain<double,4>  MChain4;

typedef ElementSet             ElemSet;
typedef ElementGroup           ElemGroup;


//=======================================================================
//   class SaturatedPorousMicroFractureModel
//=======================================================================


class SaturatedPorousMicroFractureModel : public Model
{
 public:

  typedef SaturatedPorousMicroFractureModel     Self;
  typedef Model                            Super;

  static const char*        DISP_NAMES[3];
  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;

  static const char*        GENERATE_NEW_SEED;
  static const char*        KEEP_OFF_DIAGS_PROP;
  static const char*        ARCLEN_MODE_PROP;
  static const char*        CONVEXIFY_PROP;

  static const char*        INTRIN_PERM_PROP;
  static const char*        FLUID_VISC_PROP;
  static const char*        SOLID_STIFF_PROP;
  static const char*        FLUID_STIFF_PROP;
  static const char*        POROSITY_PROP;
  static const char*        BIOT_COEFF_PROP;
  static const char*        DTIME_PROP;

  static const char*        SOLID_DENSITY_PROP;
  static const char*        FLUID_DENSITY_PROP;
  static const char*        GRAVITY_PROP;

  static const char*        FRACTURE_POROSITY_PROP;
  static const char*        RANDOM_POROSITY_PROP;
  static const char*        PERM_TYPE_PROP;

  static const char*        FIXED_ISO_PERM;
  static const char*        VOL_STRAIN_PERM;
  static const char*        TENSILE_STRAIN_PERM;
  static const char*        CUBIC_PERM;
  static const char*        CUBIC_ISO_PERM;

  static const char*        FRACTURE_TYPE_PROP;
  static const char*        GRIFFITH_ENERGY_PROP;
  static const char*        LENGTH_SCALE_PROP;
  static const char*        TENSILE_STRENGTH_PROP;
  static const char*        PENALTY_PROP;
  static const char*        PRESSURE_PSI_PROP;
  static const char*        MAX_OITER_PROP;
  static const char*        OITER_TOL_PROP;

  static const char*        BRITTLE_AT1;
  static const char*        BRITTLE_AT2;
  static const char*        LINEAR_CZM;
  static const char*        EXPONENTIAL_CZM;
  static const char*        CORNELISSEN_CZM;
  static const char*        HYPERBOLIC_CZM;

  static vector<String>     fractureModels;
  static vector<String>     permeabilityModels;

                            SaturatedPorousMicroFractureModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )  const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~SaturatedPorousMicroFractureModel  ();


 private:

  void                      getIntForce_

    ( const Vector&           force,
      const Vector&           state,
      const Vector&           state0,
      const Vector&           state00 );

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           state,
      const Vector&           state0 );

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder );

  void                      getArcFunc_

    ( const Properties&       params,
      const Properties&       globdat ); 

  void                      getUnitLoad_

    ( const Properties&       params,
      const Properties&       globdat );       

  void                      setTaylorHoodConstraints_

    ( Constraints&            cons );

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  /**
   * @brief      This function is called by getTable_.
   *
   * @param[out] table     The table
   * @param[in]  weights   The weights
   * @param[in]  contents  Filter for which data to be written to table
   * @param[in]  state     The nodal value vector
   */

  void                      getOutputData_

    ( Ref<XTable>             table,
      const Vector&           weights,
      const String&           contents,
      const Vector&           state );         

  void                      getStress_

    ( XTable&                 table,
      const Vector&           weights );

  void                      getStrain_

    ( XTable&                 table,
      const Vector&           weights );

  void                      getLocalPhaseField_

    ( XTable&                 table,
      const Vector&           weights );       

  void                      getHistory_

    ( XTable&                 table,
      const Vector&           weights );

    
  /* params.set ( "accept", false ) => keep same load and resolve
   * params.set ( "accept", true )  => advance to new load step (normal case)
   */

  void                      checkCommit_

    ( const Properties&       params,
      const Properties&       globdat );

  /* Initialize the mapping between integration points
   * and material points 
   */

  void                      initializeIPMPMap_ ();


  void                      setStepSize_

    ( const Properties&       params );

 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;

  Ref<IShape>               shape_;

  Ref<XDofSpace>            dofs_;
  IdxVector                 dofTypes_;
  IdxVector                 phiTypes_;
  IdxVector                 presTypes_;
  IdxVector                 dispTypes_;

  Matrix                    strain_;

  Ref<Material>              material_;
  Ref<PhaseFractureMaterial> phaseMaterial_;

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;

  /* 
   *  Mapping between integration points and material points  
   */

  SparseArray <int, 2>      ipMpMap_;

  /** 
   * Integer vector holds 1 or zero to indicate if element ID is 
   * active or not. Used to remove fully damaged element.
   */

  IdxVector                 isActive_;


  /*
   * Poromechanics parameters: intrinsic permeability, dynamic viscosity,
   * solid grain stiffness, fluid stiffness, porosity, biot coefficient,
   * time-step, solid density, fluid density, permeability and 
   * stabilization types, fracture porosity 
   */

  double kappa_;
  double mu_;
  double Ks_;
  double Kf_;
  double n0_;
  double alpha_;
  double dtime_;
  double rhoS_;
  double rhoF_;

  String permType_;
  bool   fracPorosity_;
  bool   randPorosity_;

  /*
   * Poromechanics derived quantities: gravitational acceleration, 
   * gravity vector, effective permeability, Voigt representation 
   * of divergence operator
   */

  bool gravity_;

  double g_;
  Vector gVec_;
  double Keff_;
  Vector voigtDiv_;

  // Derived

  Vector nrand_;

  /* Flag to off-diagonal components of the 
   * stiffness matrix is assembled or set to zero. 
   */

  bool keepOffDiags_;

  /* Flag to indicate whether this model contributes the Arc-length
   * Function
   */

  bool arcLenMode_;

  /*
   * Phase-field fracture: Fracture type, Griffith fracture energy,
   * length-scale, tensile strength, penalty 
   */

  String fracType_;
  double gc_;
  double l0_;
  double sigmaT_;
  double beta_;
  bool   pressurePsi_;
  bool   convexify_;

  /*
   * Phase-field fracture derived quantities: constants cw, p, eta,
   * a1, a2, a3, Psi0
   */

  double cw_;
  double p_;
  double eta_;
  double a1_;
  double a2_;
  double a3_;
  double Psi0_;

  double diss_;

  /*
   * Current and old stepsize (required for extrapolation)
   */

  double dt_ ;
  double dt0_;

  // Store current and old step local phase-field variable

  struct                  hist_
  {
    Flex<double>            phasef_ ;     // phase-field variable
  };

  hist_                   preHist_;      // history of the previous load step
  hist_                   newHist_;      // history of the current iteration

  // Extrapolated state vector (used for phase-field)

  Vector stateExt_;

  double errExt_;
  bool   extFail_;
  idx_t  oIter_;

  idx_t  maxOIter_;
  double oIterTol_;

};