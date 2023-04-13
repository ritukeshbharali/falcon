/** @file DirichletModel.h
 *  @brief Implements Dirichlet boundary conditions.
 *  
 *  This class implements a model for setting Dirichlet 
 *  boundary conditions for certain group of nodes. It
 *  must be included after the FE Model in *.pro file.
 *  The step size set in this model may be modified by
 *  any module that throws the action 
 *  SolverNames::SET_STEP_SIZE.
 *
 *  Author: F.P. van der Meer, f.p.vandermeer@tudelft.nl
 *          R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: May 2010
 *  
 *  Usage:      
 * 
 *    model =
 *       {
 *        type       = "Dirichlet";
 *        nodeGroups = "TopNodes";
 *        dofs       = "dy";
 *        factors    = [ -1.];
 *        stepSize   = 0.001;
 *       };
 *
 *  Updates (when, what and who)
 *     - [28 March 2022] 
 */

 /* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>

using namespace jem;

using jem::util::Properties;

using jive::Vector;
using jive::IdxVector;
using jive::StringVector;
using jive::fem::NodeSet;
using jive::fem::NodeGroup;
using jive::model::Model;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;


//=======================================================================
//   class DirichletModel
//=======================================================================

/** @brief 
 *  The DirichletModel class enforces Dirichlet boundary conditions on 
 *  user-defined group of nodes.
 */ 

class DirichletModel : public Model
{
 public:

  typedef Model             Super;
  typedef DirichletModel    Self;

  static const char*        TYPE_NAME;
  static const char*        NODE_GROUP_PROP;
  static const char*        DOFS_PROP;
  static const char*        FACTORS_PROP;
  static const char*        STEP_SIZE_PROP;

                            DirichletModel ();

                            DirichletModel

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

  virtual                  ~DirichletModel ();

 private:

  void                      init_

    ( const Properties&       globdat );

  void                      advance_

    ( const Properties&       globdat );  

  void                      applyConstraints_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      commit_

    ( const Properties&       params,
      const Properties&       globdat );
      
  void                      setStepSize_

    ( const Properties&       params );      


 private:

  StringVector              nodeGroups_;
  StringVector              dofTypes_;
  Vector                    factors_;

  Assignable<NodeSet>       nodes_;
  Ref<XDofSpace>            dofs_;
  Ref<Constraints>          cons_;

  int                       ngroups_;

  double                    stepSize_;
  double                    total0_;
  double                    total_;

};