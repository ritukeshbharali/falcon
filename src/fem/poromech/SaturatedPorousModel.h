
/** @file SaturatedPorousModel.h
 *  @brief Saturated porous media model.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 04 March 2022
 * 
 *  @note Stiffness matrix is unsymmetric, choose the linear
 *        solver accordingly.
 *
 *  Updates (when, what and who)
 *     - [11 March 2022] added option for perturbation 
 *       based numerical stabilization (RB), 
 *       DOI: 10.1002/nme2295
 * 
 *     - [14 March 2022] modified pressure equation to
 *       obtain a symmetric stiffness matrix (RB)
 * 
 *     - [22 March 2022] time step-size (dtime) can be
 *       modified by the solver that throws an action
 *       XImpActParams::SET_STEP_SIZE. An example would
 *       be the AdaptiveSteppingModule. (RB)
 * 
 *     - [05 July 2022] added getIntForce_ to compute
 *       internal force, required by quasi-Newton/line
 *       search algorithms. (RB)
 * 
 *     - [27 December 2023] removed getIntForce_,
 *       getMatrix_ returns the internal force if
 *       mbuilder = nullptr. Eliminates duplicate code. (RB)
 *       
 */

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/sparse/utilities.h>
#include <jem/util/Properties.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Globdat.h>

/* Include falcon headers */

#include "util/BasicUtils.h"
#include "materials/Material.h"

using namespace jem;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::SparseArray;
using jem::util::Properties;
using jive::Vector;
using jive::IdxVector;
using jive::IntMatrix;
using jive::Matrix;
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
//   class SaturatedPorousModel
//=======================================================================

/** @brief 
 *  The SaturatedPorousModel class implements a Finite Element Model for
 *  saturated porous media.
 * 
 *  This class implements the fully saturated porous media 
 *  model with two phases - solid and fluid. The effect of
 *  gravity is absent. The fluid pressure equation is 
 *  stabilized via perturbation (see DOI: 10.1002/nme2295).
 * 
 */

class SaturatedPorousModel : public Model
{
 public:

  typedef SaturatedPorousModel     Self;
  typedef Model                    Super;

  static const char*        DISP_NAMES[3];
  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;

  static const char*        INTRIN_PERM_PROP;
  static const char*        FLUID_VISC_PROP;
  static const char*        SOLID_STIFF_PROP;
  static const char*        FLUID_STIFF_PROP;
  static const char*        POROSITY_PROP;
  static const char*        BIOT_COEFF_PROP;
  static const char*        DTIME_PROP;

                            SaturatedPorousModel

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

  virtual                  ~SaturatedPorousModel  ();


 private:

  void                      getMatrix_

    ( Ref<MatrixBuilder>      mbuilder,
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

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;

  /* 
   *  Mapping between integration points and material points  
   */

  SparseArray <int, 2>      ipMpMap_;

  /* 
   * Integer vector holds 1 or zero to indicate if element ID is 
   * active or not. Used to remove fully damaged element.
   */

  IdxVector                 isActive_;


  /*
   * Poromechanics parameters: intrinsic permeability, dynamic viscosity,
   * solid grain stiffness, fluid stiffness, porosity, biot coefficient,
   * time-step 
   */

  double kappa_;
  double mu_;
  double Ks_;
  double Kf_;
  double phi_;
  double alpha_;
  double dtime_;

  /*
   * Poromechanics derived quantities: storage coefficient, effective
   * permeability, Voigt representation of divergence operator
   */

  double Sto_;
  double Keff_;
  Vector voigtDiv_;

};