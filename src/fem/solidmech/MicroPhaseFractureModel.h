
/** @file MicroPhaseFractureModel.h
 *  @brief Micromorphic phase-field fracture.
 *
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 14 June 2022  
 *
 *  Updates (when, what and who)
 *     - [06 August 2022] replaced hard-coded Amor phase
 *       field model with generic material update to
 *       different phase-field material models. (RB)
 * 
 *     - [25 December 2023] removed getIntForce_,
 *       getMatrix_ returns the internal force if
 *       mbuilder = nullptr. Eliminates duplicate code. (RB)
 * 
 *     - [21 November 2024] added 'keepOffDiags' option
 *       in assembling the stiffness matrix. By default,
 *       it is set to 'false'. (RB)
 */

/* Include c++ headers */

#include <vector>

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/sparse/utilities.h>
#include <jem/util/Flex.h>
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
#include <jive/implict/ArclenActions.h>

/* Include falcon headers */

#include "util/ShapeUtils.h"
#include "materials/Material.h"
#include "materials/PhaseFractureMaterial.h"

using namespace jem;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::Flex;
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
using jive::implict::ArclenActions;
using jive::implict::ArclenParams;

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
//   class MicroPhaseFractureModel
//=======================================================================

/** @brief 
 *  The MicroPhaseFractureModel class implements a micromorphic phase
 *  -field fracture FE Model (without any extrapolation technique).
 *  <a href="https://link.springer.com/article/10.1007/s00466-023-02380-1" target="_blank">Link to Article</a>
 */

class MicroPhaseFractureModel : public Model
{
  
 public:

  typedef MicroPhaseFractureModel  Self;
  typedef Model                    Super;

  static const char*        DISP_NAMES[3];
  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        RHO_PROP;
  static const char*        KEEP_OFF_DIAGS_PROP;
  static const char*        ARCLEN_MODE_PROP;

  static const char*        FRACTURE_TYPE_PROP;
  static const char*        GRIFFITH_ENERGY_PROP;
  static const char*        LENGTH_SCALE_PROP;
  static const char*        TENSILE_STRENGTH_PROP;
  static const char*        PENALTY_PROP;

  static const char*        BRITTLE_AT1;
  static const char*        BRITTLE_AT2;
  static const char*        LINEAR_CZM;
  static const char*        EXPONENTIAL_CZM;
  static const char*        CORNELISSEN_CZM;
  static const char*        HYPERBOLIC_CZM;

  static vector<String>     fractureModels;

                            MicroPhaseFractureModel

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

  virtual                  ~MicroPhaseFractureModel  ();

 private:

  void                      getMatrix_

    ( Ref<MatrixBuilder>      mbuilder,
      const Vector&           force,
      const Vector&           state );

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder );

  void                      getArcFunc_

    ( const Properties&       params,
      const Properties&       globdat );  

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

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

  void                      checkCommit_

    ( const Properties&       params );  

  /* Initialize the mapping between integration points
   * and material points 
   */

  void                      initializeIPMPMap_ ();
    
 private:

  Assignable<ElemGroup>      egroup_;
  Assignable<ElemSet>        elems_;
  Assignable<NodeSet>        nodes_;

  int                        rank_;

  Ref<IShape>                shape_;

  Ref<XDofSpace>             dofs_;
  IdxVector                  dofTypes_;
  IdxVector                  dispTypes_;
  IdxVector                  phiTypes_;

  Matrix                     strain_;

  double                     rho_;
  Ref<Material>              material_;
  Ref<PhaseFractureMaterial> phaseMaterial_;

  ShapeGradsFunc             getShapeGrads_;
  ShapeFunc                  getShapeFuncs_;

  /* 
   *  Mapping between integration points and material points  
   */

  SparseArray <int, 2>      ipMpMap_;

  /* 
   * Integer vector holds 1 or zero to indicate if element ID is 
   * active or not. Used to remove fully damaged element.
   */

  IdxVector                 isActive_;

  /* Flag to indicate off-diagonal components of the 
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

  // Store current and old step local phase-field variable

  struct                  hist_
  {
    Flex<double>          phasef_ ;     // phase-field variable
  };

  hist_                   preHist_;     // history of the previous load step
  hist_                   newHist_;     // history of the current iteration

};