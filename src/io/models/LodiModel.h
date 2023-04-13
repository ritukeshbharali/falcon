//=======================================================================
//
// Model that prints averaged load-displacement data for node group or
// nodes to file
// FPM, 29-1-2008
// Modified:
//  VP Nguyen: add evalMasterDofs(fint,cons_) so that when the nodes on which
//  external forces are applied are master nodes, their internal forces are 
//  correctly defined.
//
// Data is computed for GET_MATRIX0 and GET_INT_VECTOR
//         and printed for CHECK_COMMIT
//
//=======================================================================

#include <jem/io/PrintWriter.h>
#include <jem/io/FileWriter.h>
#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/mp/utilities.h>

#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/Assignable.h>
#include <jive/mp/Globdat.h>
#include <jive/mp/ItemMask.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/fem/ElementSet.h>
#include <jive/fem/NodeGroup.h>

using namespace jem;
using jem::io::endl;

using jem::util::Properties;
using jem::io::PrintWriter;
using jem::io::FileWriter;
using jem::mp::allsum;
using jem::mp::Context;
using jive::Vector;
using jive::IdxVector;
using jive::IdxMatrix;
using jive::StringVector;
using jive::util::DofSpace;
using jive::util::Assignable;
using jive::util::evalMasterDofs;
using jive::util::Constraints;
using jive::model::Model;
using jive::fem::NodeSet;
using jive::fem::NodeGroup;
using jive::fem::ElementSet;


//=======================================================================
//   class LodiModel
//=======================================================================

/** @brief 
 *  The LodiModel class computes the load and displacement on a
 *  user-defined group of nodes.
 */ 

class LodiModel : public Model
{
 public:

  static const char*        TYPE_NAME;
  static const char*        NODES;
  static const char*        GROUP;
  static const char*        FILE_NAME;


                            LodiModel

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~LodiModel  ();


 private:

  Ref<Context>              mpx_;
  Ref<DofSpace>             dofs_;
  Ref<Constraints>          cons_;
  Ref<PrintWriter>          out_;
  Assignable<ElementSet>    elems_;
  Assignable<NodeSet>       nodes_;

  IdxVector                 inodes_;
  IdxMatrix                 idofs_;

  idx_t                     nn_;
  idx_t                     ndof_;

  String                    file_;

  bool                      toFile_;
  bool                      evalMasterDof_;
};


