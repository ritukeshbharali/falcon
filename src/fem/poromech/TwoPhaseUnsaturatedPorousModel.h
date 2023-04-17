
/** @file TwoPhaseUnsaturatedPorousModel.h
 *  @brief Implements unsaturated two-phase porous media model.
 *  
 *  This class implements a mass conserving unsaturated two-phase
 *  porous media model with two phases - solid and fluid. The 
 *  fluid pressure equation is stabilized via perturbation 
 *  (see DOI: 10.1002/nme2295). The gas is assumed to be at 
 *  atmospheric pressure, equal to zero.
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
 */

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
#include "materials/RetentionMaterial.h"

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
//   class TwoPhaseUnsaturatedPorousModel
//=======================================================================


class TwoPhaseUnsaturatedPorousModel : public Model
{
 public:

  typedef TwoPhaseUnsaturatedPorousModel     Self;
  typedef Model                              Super;

  static const char*        DISP_NAMES[3];
  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        RETENTION_PROP;

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

                            TwoPhaseUnsaturatedPorousModel

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

  virtual                  ~TwoPhaseUnsaturatedPorousModel  ();


 private:

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           state,
      const Vector&           state0 );

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder );

  void                      setConstraints_

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

  void                      getHistory_

    ( XTable&                 table,
      const Vector&           weights );

    
  /* params.set ( "accept", false ) => keep same load and resolve
   * params.set ( "accept", true )  => advance to new load step (normal case)
   */

  void                      checkCommit_

    ( const Properties&       params );

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
  IdxVector                 presTypes_;
  IdxVector                 dispTypes_;

  Matrix                    strain_;

  Ref<Material>             material_;
  Ref<RetentionMaterial>    retention_;

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


  /**
   * Poromechanics parameters: intrinsic permeability, dynamic viscosity,
   * solid grain stiffness, fluid stiffness, porosity, biot coefficient,
   * time-step, solid density, fluid density 
   */

  double kappa_;
  double mu_;
  double Ks_;
  double Kf_;
  double phi_;
  double alpha_;
  double dtime_;
  double rhoS_;
  double rhoF_;

  /**
   * Poromechanics derived quantities: gravitational acceleration, 
   * gravity vector, effective permeability, Voigt representation 
   * of divergence operator
   */

  bool gravity_;

  double g_;
  Vector gVec_;
  double Keff_;
  Vector voigtDiv_;

};