
/** @file LinearElasticityModel.cpp
 *  @brief Implements the linear elasticity model.
 *  
 *  This class implements a finite element model with
 *  linear elastic material law.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 02 March 2022
 *
 *  Updates (when, what and who)
 *     - [25 March 2022] added option for integration
 *       point averaged element data. getElemStress_
 *       and getElemStrain_ are the new functions.
 *       (RB)
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
//   class LinearElasticityModel
//=======================================================================

/** @brief 
 *  The LinearElasticityModel class implements a Finite Element Model to
 *  be used with linear elastic materials.
 */

class LinearElasticityModel : public Model
{
  
 public:

  typedef LinearElasticityModel  Self;
  typedef Model                  Super;

  static const char*        DOF_NAMES[3];
  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;
  static const char*        RHO_PROP;

                            LinearElasticityModel

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

  virtual                  ~LinearElasticityModel  ();

 private:

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           state );

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder );

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

  void                      getElemStress_

    ( XTable&                 table,
      const Vector&           weights );  

  void                      getStrain_

    ( XTable&                 table,
      const Vector&           weights );

  void                      getElemStrain_

    ( XTable&                 table,
      const Vector&           weights );      

  void                      getHistory_

    ( XTable&                 table,
      const Vector&           weights );

  void                      checkCommit_

    ( const Properties&       params );

  /* Initialize the mapping between integration points
   * and material points 
   */

  void                      initializeIPMPMap_ ();


 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;

  Ref<IShape>               shape_;

  Ref<XDofSpace>            dofs_;
  IdxVector                 dofTypes_;

  Matrix                    strain_;

  double                    rho_;
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

};