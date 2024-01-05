
/** @file ProstateCancerModel.h
 *  @brief Implements a prostate cancer model.
 *  
 *  This class implements a phase-field model for prostate
 *  case (see DOI: 10.1073/pnas.1615791113).
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 24 November 2022
 * 
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022] 
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

/* Include user-defined class headers */

#include "util/ShapeUtils.h"

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
//   class ProstateCancerModel
//=======================================================================


class ProstateCancerModel : public Model
{
 public:

  typedef ProstateCancerModel     Self;
  typedef Model                   Super;

  static const char*        SHAPE_PROP;
  static const char*        DTIME_PROP;
  static const char*        LAMBDA_PROP;
  static const char*        TAU_PROP;
  static const char*        GROWTH_RATE_PROP;
  static const char*        APOPTOSIS_RATE_PROP;
  static const char*        EPSILON_PROP;
  static const char*        NUTRIENT_AVG_PROP;
  static const char*        NUTRIENT_DEV_PROP;
  static const char*        NUTRIENT_CONSUMPTION_PROP;
  static const char*        NUTRIENT_DECAY_PROP;
  static const char*        GENERATE_NEW_SEED;
  
                            ProstateCancerModel

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

  virtual                  ~ProstateCancerModel  ();


 private:

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           state,
      const Vector&           state0 );


  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getElemSource_

    ( XTable&                 table,
      const Vector&           weights );


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
  IdxVector                 sigTypes_;

  ShapeGradsFunc            getShapeGrads_;
  ShapeFunc                 getShapeFuncs_;

  /* 
   *  Mapping between integration points and material points  
   */

  SparseArray <int, 2>      ipMpMap_;

  // Input parameters

  double dtime_;
  double lambda_;
  double tau_;
  double chi_;
  double A_;
  double eps_;
  double delta_;
  double gamma_;

  double savg_;
  double sdev_;

  // Derived

  Vector srand_;

};


