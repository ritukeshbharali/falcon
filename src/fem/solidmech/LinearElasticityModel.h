
/** @file LinearElasticityModel.h
 *  @brief Linear elasticity model.
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
 *  The LinearElasticityModel class implements a linear elasticity FE 
 *  model to be used with linear elastic materials (Hooke, Orthotropic).
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

  /**
   * @brief getMatrix_ function assembles the stiffness matrix 
   *        and the internal force vector.
   *
   * @param[out] mbuilder  Ref to Matrix Builder (Stiffness Matrix)
   * @param[in]  force     Internal force vector
   * @param[in]  state     Current state (solution) vector
   */

  void                      getMatrix_

    ( MatrixBuilder&          mbuilder,
      const Vector&           force,
      const Vector&           state );

  /**
   * @brief getMatrix2_ function assembles the mass matrix.
   *
   * @param[out] mbuilder  Ref to Matrix Builder (Mass Matrix)
   */

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder );

  /**
   * @brief getTable_ function returns post-processing quantities
   *        in a tabular format.
   *
   * @param[in,out] params    Properties object containing the table
   * @param[in]  globdat      Global database
   */  

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  /**
   * @brief getOutputData_ function returns post-processing quantities
   *        in a tabular format (originally written by Frans van der Meer).
   *
   * @param[out] table     Table containing post-processing quantities
   * @param[out] weights   Weights
   * @param[out] contents  Specific contents requested by the user.
   * @param[in]  state     Current state (solution) vector
   */  

  void                      getOutputData_

    ( Ref<XTable>             table,
      const Vector&           weights,
      const String&           contents,
      const Vector&           state );

  /**
   * @brief getStress_ function returns nodal stresses in a table along 
   *        with the weights.
   *
   * @param[out] table     Table containing post-processing quantities
   * @param[out] weights   Weights
   */

  void                      getStress_

    ( XTable&                 table,
      const Vector&           weights );

  /**
   * @brief getElemStress_ function returns element (average of 
   *        integration points) stresses in a table along with the weights.
   *
   * @param[out] table     Table containing post-processing quantities
   * @param[out] weights   Weights
   */

  void                      getElemStress_

    ( XTable&                 table,
      const Vector&           weights );

  /**
   * @brief getStrain_ function returns nodal strains in a table along 
   *        with the weights.
   *
   * @param[out] table     Table containing post-processing quantities
   * @param[out] weights   Weights
   */     

  void                      getStrain_

    ( XTable&                 table,
      const Vector&           weights );
  
  /**
   * @brief getElemStrain_ function returns element (average of
   *        integration points) strains in a table along with the weights.
   *
   * @param[out] table     Table containing post-processing quantities
   * @param[out] weights   Weights
   */

  void                      getElemStrain_

    ( XTable&                 table,
      const Vector&           weights );

  /**
   * @brief getHistory_ function does nothing for linear elastic models.
   *
   * @param[out] table     Table containing post-processing quantities
   * @param[out] weights   Weights
   */  

  void                      getHistory_

    ( XTable&                 table,
      const Vector&           weights );

  /**
   * @brief checkCommit_ function does nothing for linear elastic models.
   *
   * @param[in,out] params Properties
   */  

  void                      checkCommit_

    ( const Properties&       params );

  /**
   * @brief initializeIPMPMap_ function initializes the mapping between 
   * integration points and material points.
   */

  void                      initializeIPMPMap_ ();


 private:

  Assignable<ElemGroup>     egroup_;         /**< Element Group assigned to this model */
  Assignable<ElemSet>       elems_;          /**< Element Set corresponding to Element Group*/
  Assignable<NodeSet>       nodes_;          /**< (Unique) Node Set corresponding to Element Set*/

  int                       rank_;           /**< Mesh rank */
  Ref<IShape>               shape_;          /**< Ref to Shape function */

  Ref<XDofSpace>            dofs_;           /**< Ref to the DoF space */
  IdxVector                 dofTypes_;       /**< DoF types vector */ 

  Matrix                    strain_;         /**< Matrix that stores strain */

  double                    rho_;            /**< Density (reqd. for mass matrix) */
  Ref<Material>             material_;       /**< Ref to the material model */

  ShapeGradsFunc            getShapeGrads_;  /**< Function for rank-based B-matrix */
  ShapeFunc                 getShapeFuncs_;  /**< Function for rank-based N-matrix */

  SparseArray <int, 2>      ipMpMap_;        /**< Mapping between integration and material points */
  IdxVector                 isActive_;       /**< Vector that removes fully damaged elements */

};