/*
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements the gradient enhanced damage model
 *  using finite elements. The number of elements and the 
 *  integration scheme/element can be changed during the 
 *  simulation.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 */

/** @file GradientDamageModel.h
 *  @brief Gradient enhanced damage model.
 * 
 *  Original Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *  
 *  Updates (when, what and who)
 *     - [21 February 2024] split original code into
 *     header and source files for readability, modified
 *     getMatrix_ to return only internal force when
 *     mbuilder = nullptr, added functions to obtain the
 *     strains, stresses, and damage. (RB)
 */


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

#include "util/ShapeUtils.h"
#include "materials/Material.h"
#include "materials/DamageMaterial.h"

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
//   class GradientDamageModel
//=======================================================================

/** @brief 
 *  The GradientDamageModel class implements the gradient enhanced damage
 *  FE Model proposed in Peerling et. al. (1996). 
 *  <a href="https://doi.org/10.1002/(SICI)1097-0207(19961015)39:19%3C3391::AID-NME7%3E3.0.CO;2-D" target="_blank">Link to Article</a>
 */

class GradientDamageModel : public Model
{
 public:

  typedef GradientDamageModel     Self;
  typedef Model                   Super;

  static const char*        DISP_NAMES[3];

  static const char*        SHAPE_PROP;
  static const char*        MATERIAL_PROP;


                            GradientDamageModel

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

  virtual                  ~GradientDamageModel  ();


 private:

  void                      getMatrix_

    ( Ref<MatrixBuilder>      mbuilder,
      const Vector&           force,
      const Vector&           state );

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

  void                      getDamage_

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

 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;

  Ref<IShape>               shape_;

  Ref<XDofSpace>            dofs_;
  IdxVector                 dofTypes_;
  IdxVector                 eqvTypes_;
  IdxVector                 dispTypes_;

  Matrix                    strain_;

  Ref<DamageMaterial>       material_;

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