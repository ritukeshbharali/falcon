
/** @file ParaviewModule.h
 *  @brief Implements module to write VTK files in parallel.
 *  
 *  This class implements a module to write VTK output 
 *  (*.vtu files) in parallel. Every process prints the
 *  mesh and data it owns. Solution vector (STATE) is
 *  printed as a default choice. User can request other
 *  data such as stress, strain etc. Two types of data
 *  may be obtained (i) nodal, and (ii) element. Nodal
 *  data performs a nodal averaging and is continuous
 *  (smoothened). Element data is averaged over the
 *  integration points (since VTK CellData allows only
 *  one value per cell).
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 25 March 2022
 * 
 *  NOTE: This model must be included after the Solver
 *  Module in the *.pro file. Data is printed every
 *  'interval' step(s).
 * 
 *  Updates (when, what and who)
 *  - [28 March 2022] throws an IllegalInput exception
 *       when dataType_ is not nodes or elems. (RB)
 *  - [18 June 2022] ignores solution fields defined on
 *    dummy nodes, the solution fields are printed per
 *    components instead of the component wise merged
 *    style. Solution fields on the dummy nodes may be
 *    extracted as PointData or CellData. (RB)
 *       
 */


#ifndef PARAVIEW_MODULE_H
#define PARAVIEW_MODULE_H

/* Include jem and jive headers */

#include <jive/app/Module.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Assignable.h>
#include <jem/mp/utilities.h>
#include <jive/mp/Globdat.h>
#include <jem/util/ArrayBuffer.h>

#include <map>


namespace jem
{
  namespace io
  {
    class PrintWriter;	  
  }
}

namespace jive
{
  namespace util
  {
    class XTable;	  
    class XDofSpace;
  }
}

using namespace jem;
using jem::String;
using jem::io::PrintWriter;
using jem::util::Properties;
using jem::mp::Context;

using jive::app::Module;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::util::Assignable;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::IdxVector;
using jive::StringVector;
using jem::util::ArrayBuffer;



//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef ElementSet              ElemSet;
typedef ElementGroup            ElemGroup;


//=======================================================================
//   class ParaviewModule
//=======================================================================

/** @brief 
 *  The \c ParaviewModule class writes vtk files for post-processing 
 *  simulation data in Paraview, MayaVi or VisIt.
 * 
 *  Below is an example how \c ParaviewModule may be used:
 * 
 *  \code
 *  vtk = "paraview"
 *    {
 *       fileName = "$(CASE_NAME)_out";
 *       elements = "DomainElems";
 *       printInterval = 1;
 *       pointData  = ["stress","strain"];
 *       cellData  = ["stress","strain"];
 *    };
 *  \endcode
 */ 

class ParaviewModule : public Module
{
 public:

  typedef ParaviewModule    Self;
  typedef Module            Super;

  static const char*        FILENAME_PROP;
  static const char*        PRINT_INTERVAL_PROP;
  static const char*        POINT_DATA_PROP;
  static const char*        CELL_DATA_PROP;

  explicit                  ParaviewModule

    ( const String&           name = "vtk" );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual Status            run

    ( const Properties&       globdat );

  virtual void              shutdown

    ( const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )        const;

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       globdat,
      const Properties&       props );

 protected:

  virtual                  ~ParaviewModule   ();

 private:

  // Write PVD file

  void                     writePVDHeader_  ( Ref<PrintWriter>  pvdWriter  );
  void                     writePVDClosure_ ( Ref<PrintWriter>  pvdWriter  );

  void                     writePVDTimeStamp_ ( Ref<PrintWriter>  pvdWriter,
                                                const String&     iFileName,
                                                const double&     iTime     );

  // Write VTU file

  void                     writeVTUHeader_  ( Ref<PrintWriter>  vtuWriter  );
  void                     writeVTUClosure_ ( Ref<PrintWriter>  vtuWriter  );

  // Write FE information

  void                     writeMesh_       ( Ref<PrintWriter>  vtuWriter  );

  void                     writeState_      ( Ref<PrintWriter>  vtuWriter,
                                              const Properties& globdat  );

  // Write tabular data

  void                     writeTable_      ( Ref<PrintWriter>  vtuWriter,
                                              Ref<XTable>       table,
                                              const String&     tname,
                                              const int         rowCount );

 private:

  Ref<Context>             mpx_;

  String                   fileName_;

  int                      rank_;

  // Mesh details

  Assignable<ElemGroup>    egroup_;
  Assignable<ElemSet>      elems_;
  Assignable<NodeSet>      nodes_;
  IdxVector                inodes_;
  IdxVector                ielems_;
  int                      elemCount_;
  int                      nodeCount_;

  // Element mapping to VTK format
  
  int                      nVtxPerElem_;
  IdxVector                vtxNodes_;
  int                      vtkCellCode_;
  int                      nodesPerElem_;

  // Dofspace details

  Ref<XDofSpace>           dofs_;
  StringVector             dofNames_;
  IdxVector                dofTypes_;
  idx_t                    dofTypeCount_;

  // Printing request details

  int                      interval_;
  StringVector             pointData_;
  StringVector             cellData_;

  ArrayBuffer<double>      timeStamps_;
  ArrayBuffer<double>      stepStamps_;

};


#endif
