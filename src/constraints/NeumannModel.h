
/** @file NeumannModel.h
 *  @brief Implements Neumann boundary conditions.
 *  
 *  This class implements the computation of traction
 *  loads on the Neumann boundary. It is compatible 
 *  with any FE model, since does not require any pre-
 *  knowledge of the dofs. The dofs are extracted from
 *  the (X)Dofspace. This model must be included after
 *  the FEmodel in the *.pro file.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 13 March 2022
 *
 *  Usage:      
 * 
 *    model =
 *       {
 *        type     = "Neumann";
 *        elements = "TopElems";
 *        loads    = [0.0,-1.0e+04,0.0]; // size = dofs
 *        shape  =
 *         {
 *           type  = "BLine2";
 *           shapeFuncs  =
 *           {
 *             type  = "Linear";
 *           };
 *           intScheme = "Gauss1";
 *         };
 *       };
 *
 *  Updates (when, what and who)
 *     - [22 March 2022] throws an IllegalInputException
 *       when size of loads is not equal to dofs.
 *       
 */

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Properties.h>
#include <jive/util/Assignable.h>
#include <jive/util/utilities.h>
#include <jive/util/XDofSpace.h>
#include <jive/geom/Line.h>
#include <jive/geom/BoundaryShape.h>
#include <jive/geom/BShapeFactory.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/fem/ElementGroup.h>
#include <jive/femodel/Names.h>

using namespace jem;

using jem::util::Properties;
using jive::Vector;
using jive::StringVector;
using jive::Matrix;
using jive::IdxVector;
using jive::util::Assignable;
using jive::util::XDofSpace;
using jive::model::Model;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::geom::BShape;


//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef ElementSet    ElemSet;
typedef ElementGroup  ElemGroup;


//=======================================================================
//   class NeumannModel
//=======================================================================

/** @brief 
 *  The NeumannModel class enforces Neumann boundary conditions on 
 *  user-defined group of elements.
 */ 

class NeumannModel : public Model
{
 public:

  typedef Model             Super;
  typedef NeumannModel      Self;

  static const char*        TYPE_NAME;
  static const char*        LOADS_PROP;
  static const char*        SHAPE_PROP;

                            NeumannModel ();

                            NeumannModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~NeumannModel ();

 private:

  void                      getExtForce_

    ( const Vector&           fext,
      double                  scale,
      const Properties&       globdat )      const;

 private:

  Assignable<ElemGroup>     egroup_;
  Assignable<ElemSet>       elems_;
  Assignable<NodeSet>       nodes_;

  int                       rank_;

  Ref<BShape>               shape_;

  Ref<XDofSpace>            dofs_;
  StringVector              dofNames_;
  IdxVector                 dofTypes_;

  Vector                    loads_;
};