
/** @file PeriodicModel.h
 *  @brief Implements Periodic boundary conditions.
 *  
 *  This class implements a model for setting Periodic 
 *  boundary conditions for certain group of nodes. It
 *  must be included after the FE Model in *.pro file.
 *  The step size set in this model may be modified by
 *  any module that throws the action 
 *  SolverNames::SET_STEP_SIZE.
 * 
 *  Equation: nodeGroup1 = nodeGroup0 + scale * stepsize
 *  Jive:       idof        jdof      + coeff * rval            
 */

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/util/Properties.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/util/DofSpace.h>
#include <jive/util/XDofSpace.h>

/* Include user-defined class headers */

#include "util/XNames.h"

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
using jive::util::DofSpace;
using jive::util::XDofSpace;


//=======================================================================
//   class PeriodicModel
//=======================================================================


class PeriodicModel : public Model
{
 public:

  typedef Model             Super;
  typedef PeriodicModel     Self;

  static const char*        TYPE_NAME;
  static const char*        NODE_GROUPS0_PROP;
  static const char*        NODE_GROUPS1_PROP;
  static const char*        DOFS_PROP;
  static const char*        FACTORS_PROP;
  static const char*        STEP_SIZE_PROP;

                            PeriodicModel ();

                            PeriodicModel

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

  virtual                  ~PeriodicModel ();

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

  StringVector              nodeGroups0_;
  StringVector              nodeGroups1_;
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