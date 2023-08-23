
/** @file ParaviewModule.cpp
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
 *  Usage:      
 * 
 *  vtk = "vtkWriter"
 *    {
 *       fileName = "$(CASE_NAME)_out";
 *       elements = "DomainElems";
 *       interval = 1;
 *       data     = ["stress","strain"];
 *       dataType = "nodes";   // "elems"
 *    };
 *
 *  Updates (when, what and who)
 *     - [28 March 2022] throws an IllegalInput exception
 *       when dataType_ is not nodes or elems.
 *       
 */


/* Include jem and jive headers */

#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/base/limits.h>
#include <jem/io/FileWriter.h>
#include <jem/io/PrintWriter.h>
#include <jem/util/Properties.h>
#include <jive/util/Globdat.h>
#include <jive/util/ItemSet.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/DenseTable.h>
#include <jive/util/Printer.h>
#include <jive/model/StateVector.h>
#include <jive/model/Model.h>
#include <jive/model/Actions.h>
#include <jive/implict/SolverInfo.h>
#include <jive/app/ModuleFactory.h>

/* Include user-defined class headers */

#include "ParaviewModule.h"
#include "util/BasicUtils.h"
#include "util/XNames.h"

#include <stdio.h>

using namespace jem;
using jem::io::FileWriter;
using jive::Matrix;
using jive::IdxVector;
using jive::Vector;
using jem::io::endl;
using jive::util::XDofSpace;
using jive::model::StateVector;
using jive::util::Globdat;
using jive::util::Printer;
using jive::util::DenseTable;


//=======================================================================
//   class ParaviewModule
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  ParaviewModule::FILENAME_PROP        = "fileName";
const char*  ParaviewModule::PRINT_INTERVAL_PROP  = "printInterval";
const char*  ParaviewModule::POINT_DATA_PROP      = "pointData";
const char*  ParaviewModule::CELL_DATA_PROP       = "cellData";



//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


ParaviewModule::ParaviewModule ( const String& name ) :
  Super ( name )
{}

ParaviewModule::~ParaviewModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status ParaviewModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
   using jive::mp::Globdat;

   // MPI ( Get the number of process )  

   mpx_    = Globdat::getMPContext ( globdat );
   nProcs_ = 1;                                   // maybe not required!
   nProcs_ = mpx_-> size ();

   Properties  myProps = props.getProps ( myName_ );
   Properties  myConf  = conf.makeProps ( myName_ );

   const String context = getContext();

   // Get the element group and extract corresponding nodes and dof space

   egroup_ = ElemGroup::get      ( myConf, myProps, globdat, context );
   elems_  = egroup_.getElements ( );
   nodes_  = elems_ .getNodes    ( );
   dofs_   = XDofSpace::get      ( nodes_.getData(), globdat );

   // Get some additional data ( number of nodes, elements, mesh rank,
   // dofs per node )
  
   nodeCount_    = nodes_.size     ( );
   elemCount_    = egroup_.size    ( );
   rank_         = nodes_.rank     ( );
   dofTypeCount_ = dofs_->typeCount( );

   // Get element information (element indices and nodes per element)

   ielems_.resize ( elemCount_ );
   ielems_        = egroup_.getIndices ();
   noNodePerElem_ = elems_.getElemNodeCount(ielems_[0]);

   // Get number of dofs from the XDofSpace

   dofNames_.resize ( dofTypeCount_ );
   dofTypes_.resize ( dofTypeCount_ );

   dofNames_ = dofs_->getTypeNames();

   for ( idx_t i = 0; i < dofTypeCount_; i++ )
   {
     dofTypes_[i] = dofs_->addType ( dofNames_[i] );
   }

   // Get element types

   // (2D)

   if ( rank_ == 2 )
   {
     // Check for Triangles (T3, T6)

     if ( noNodePerElem_ == 3 )
     {
       numVertexesPerCell_ = 3;
       vtkCellCode_        = 5;
     }
     else if ( noNodePerElem_ == 6 )
     {
       numVertexesPerCell_ = 3;
       vtkCellCode_        = 5;

       // Get corner/vertex nodes for T6
       cornerNodes_.resize( numVertexesPerCell_ );
       cornerNodes_ = {0,2,4};
     }

     // Check for Quads (Q4, Q8, Q9)

     else if ( noNodePerElem_ == 4 )
     {
       numVertexesPerCell_ = 4;
       vtkCellCode_        = 9;
     }
     else if ( noNodePerElem_ == 8 || noNodePerElem_ == 9 )
     {
       numVertexesPerCell_ = 4;
       vtkCellCode_        = 9;

       // Get corner/vertex nodes for Q8, Q9
       cornerNodes_.resize( numVertexesPerCell_ );
       cornerNodes_ = {0,2,4,6};
     }

     else
     {
       const String  context = getContext ();

       throw IllegalInputException (
         context,
         String::format (
           "invalid element type: (should be T3, T6, Q4, Q8 or Q9)"
        )
       );
     }
   }

   // (3D)

   else if ( rank_ == 3 )
   {
     // Check for Tetrahedrons (Tet4, Tet10)

     if ( noNodePerElem_ == 4 )
     {
       numVertexesPerCell_ = 4;
       vtkCellCode_        = 10;
     }
     else if ( noNodePerElem_ == 10 )
     {
       numVertexesPerCell_ = 4;
       vtkCellCode_        = 10;
       
       // Get corner/vertex nodes for Tet10
       cornerNodes_.resize( numVertexesPerCell_ );
       cornerNodes_ = {0,2,4,9}; 
     }

     // Check for Hexahedrons (Hex8, Hex20)

     else if ( noNodePerElem_ == 8 )
     {
       numVertexesPerCell_ = 8;
       vtkCellCode_        = 12;
     }
     else if ( noNodePerElem_ == 20 )
     {
       numVertexesPerCell_ = 8;
       vtkCellCode_        = 12;

       // Get corner/vertex nodes for Hex20
       cornerNodes_.resize( numVertexesPerCell_ );
       cornerNodes_ = {0,2,4,6,12,14,16,18};
     }

     else
     {
       const String  context = getContext ();

       throw IllegalInputException (
         context,
         String::format (
           "invalid element type: (should be Tet4, Tet10, Hex8 or Hex20)"
        )
       );
     }

   }

   return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------

Module::Status ParaviewModule::run ( const Properties& globdat )
{
  using jem::Ref;
  using jem::System;
  using jive::Vector;
  using jive::util::ItemSet;
  using jive::model::Model;
  using jive::model::ActionParams;
  using jive::model::Actions;

  // Get context

  const String context = getContext ();

  // Get current MPI rank ( every process prints the mesh/data it owns )

  int myRank = mpx_-> myRank ();

  // Get the (converged) step. At this point, the solver has already
  // updated step by +1. This change (+1) is rolled back for printing
  // output files ( pvd, vtu ).

  int step; globdat.get ( step, Globdat::OLD_TIME_STEP ); step -= 1;

  // Post-processing data is printed every 'interval_' steps
  // By default, data from the first step, i.e, '0' is always 
  // printed.

  if ( ( step % interval_) != 0 && step > 0 ) return OK;

  // Post-processing is only done when the solution is accepted 
  
  bool accepted = true;
  globdat.find ( accepted, "var.accepted" );

  if ( !accepted ) return OK;

  // Get the current time ( for time-dependent problems, 
  // time stamps are stored in pvd file )

  double time;

  if ( globdat.find ( time, Globdat::OLD_TIME ) )
  {
    timeStamps_.pushBack ( time );
  }
  else
  {
    timeStamps_.pushBack ( step );
  }

  stepStamps_.pushBack ( step );

  // -----------------------------------------------------------------------------
  // Print a PVD file ( overwritten on every step )
  // -----------------------------------------------------------------------------

  // Allocate string to store pvd filename

  const String fileNamePVD = fileName_ + ".pvd";

  // Initialize a PVD writer

  Ref<PrintWriter> pvdWriter = newInstance<PrintWriter> ( 
                               newInstance<FileWriter>  ( fileNamePVD ) );

  // Write PVD header

  writePVDHeader_ ( pvdWriter );

  // Print timestamps

  for ( idx_t iStep = 0; iStep < timeStamps_.size(); iStep++ )
  {
    const String iFileNameVTU = fileName_ + "_p" + String ( myRank )
                                          + "_"  + String ( (int) stepStamps_[iStep]  )
                                          + ".vtu";

    writePVDTimeStamp_ ( pvdWriter, iFileNameVTU, timeStamps_[iStep] );
  }

  // Write PVD closure, flush and close the file

  writePVDClosure_ ( pvdWriter );

  pvdWriter->flush();
  pvdWriter->close();

  // -----------------------------------------------------------------------------
  // Print a VTU file ( for current step )
  // -----------------------------------------------------------------------------

  const String fileNameVTU = fileName_ + "_p" + String ( myRank )
                                       + "_"  + String ( step   )
                                       + ".vtu";

  // Initialize a VTU writer

  Ref<PrintWriter> vtuWriter = newInstance<PrintWriter> ( 
                               newInstance<FileWriter>  ( fileNameVTU ) );

  // Write VTU header and mesh

  writeVTUHeader_ ( vtuWriter );
  writeMesh_      ( vtuWriter );

  // Begin printing PointData

  *vtuWriter << "<PointData  Vectors=\"NodalData\"> \n";

  // Print the solution ( in Jive terminology, 'state' )

  writeState_ ( vtuWriter, globdat );

  // Check if additional PointData is requested

  if ( pointData_.size() > 0 )
  {
    Properties      params  ("actionParams");
    Vector          weights ( nodeCount_ );

    for ( idx_t iData = 0; iData < pointData_.size(); iData++ )
    {
      String               tname ( pointData_[iData] );
      weights            = 0.0;

      Ref<Model>  model  = Model::get ( globdat, getContext() );
      Ref<XTable> xtable = newInstance<DenseTable>  ( tname, nodes_.getData() );

      params.set ( ActionParams::TABLE,         xtable  );
      params.set ( ActionParams::TABLE_NAME,    tname   );
      params.set ( ActionParams::TABLE_WEIGHTS, weights );

      bool updated = model->takeAction ( Actions::GET_TABLE, params, globdat );
      
      weights = where ( abs( weights ) < Limits<double>::TINY_VALUE, 1.0, 1.0 / weights );
      xtable -> scaleRows    ( weights );

      if ( updated ) 
      {
        writeTable_   ( vtuWriter, xtable, tname );
      }
    }
  }

  // Close PointData

  *vtuWriter << "</PointData> \n";

  // Check if additional CellData is requested

  if ( cellData_.size() > 0 )
  {
    Properties      params  ("actionParams");
    Vector          weights ( elemCount_   );  // probably not required!

    // Begin printing CellData

    *vtuWriter << "<CellData> \n";

    for ( idx_t iData = 0; iData < cellData_.size(); iData++ )
    {
      String               tname ( cellData_[iData] );
      weights            = 1.0;

      Ref<Model>  model  = Model::get ( globdat, getContext() );
      Ref<XTable> xtable = newInstance<DenseTable>  ( tname, elems_.getData() );

      params.set ( ActionParams::TABLE,         xtable  );
      params.set ( ActionParams::TABLE_NAME,    tname   );
      params.set ( ActionParams::TABLE_WEIGHTS, weights );

      bool updated = model->takeAction ( Actions::GET_TABLE, params, globdat );

      if ( updated ) 
      {
        writeTable_   ( vtuWriter, xtable, tname, elemCount_ );
      }
    }

    // Close CellData

    *vtuWriter << "</CellData> \n";

  }

  // Write Closure

  writeVTUClosure_   ( vtuWriter );

  // Flush and Close file

  vtuWriter->flush();
  vtuWriter->close();

  return OK;

}


//-----------------------------------------------------------------------
//   writePVDHeader_
//-----------------------------------------------------------------------
    

void ParaviewModule::writePVDHeader_ ( Ref<PrintWriter> pvdWriter )
{
  *pvdWriter << "<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n";
  *pvdWriter << "<Collection>\n";
}


//-----------------------------------------------------------------------
//   writePVDClosure_
//-----------------------------------------------------------------------
    

void ParaviewModule::writePVDClosure_ ( Ref<PrintWriter> pvdWriter )
{
  *pvdWriter << "</Collection>\n";
  *pvdWriter << "</VTKFile>\n";
}


//-----------------------------------------------------------------------
//   writePVDTimeStamp_
//-----------------------------------------------------------------------
    

void ParaviewModule::writePVDTimeStamp_ ( Ref<PrintWriter> pvdWriter,
                                          const String&    iFileName,
                                          const double&    iTime     )
{
  *pvdWriter << "<DataSet file='" + iFileName 
               + "' groups='' part='0' timestep='"
               + String(iTime)+"'/>\n";
}


//-----------------------------------------------------------------------
//   writeVTUHeader_
//-----------------------------------------------------------------------
    

void ParaviewModule::writeVTUHeader_ ( Ref<PrintWriter> vtuWriter )
{
  const String str1 = String::format("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"> \n"
                                 "<UnstructuredGrid> \n"
                                   "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\"> \n", 
                             nodeCount_, elemCount_ );

  *vtuWriter << str1;
}



//-----------------------------------------------------------------------
//   writeMesh_
//-----------------------------------------------------------------------
    

void ParaviewModule::writeMesh_ ( Ref<PrintWriter> vtuWriter )
{

  // Write node coordinates

  *vtuWriter << "<Points> \n";
  *vtuWriter << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" >\n";

  Matrix coords(rank_,nodeCount_);

  nodes_.getCoords ( coords );

  if ( rank_ == 2 )
  {
    for(idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuWriter << coords(0,i) << " " << coords(1,i) << " " << 0. << '\n';
    }
  }
  else
  {
    for(idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuWriter << coords(0,i) << " " << coords(1,i) << " " << coords(2,i) << '\n';
    }
  }

  *vtuWriter << "</DataArray>\n </Points>\n";

  // Write element connectivity

  IdxVector inodes(noNodePerElem_);

  *vtuWriter << "<Cells> \n";
  *vtuWriter << "<DataArray  type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  { 
    int  ielem = ielems_[ie]; 
    elems_.getElemNodes  ( inodes, ielem );

    // Write only the corner nodes for higher order elements
    // Middle nodes are not required for mesh visualization

    if (noNodePerElem_ == numVertexesPerCell_ )
    {
      for ( idx_t in = 0; in < noNodePerElem_; in++ ){
        *vtuWriter << inodes[in] << " "; 
      }
    }
    else
    {
      for ( idx_t in = 0; in < numVertexesPerCell_; in++ ){
        *vtuWriter << inodes[cornerNodes_[in]] << " "; 
      }
    }

    *vtuWriter << endl;
  }

  *vtuWriter << "</DataArray>\n";

  // write cell offset

  *vtuWriter << "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\"> \n";

  int offset = 0;
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    offset += numVertexesPerCell_;
    *vtuWriter <<  offset << endl;
  }

  *vtuWriter << "</DataArray>\n";

  // Print cell types
  
  *vtuWriter << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\"> \n";
  
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
     *vtuWriter <<  vtkCellCode_ << endl;
  }

  *vtuWriter << "</DataArray> \n </Cells> \n";

}

//-----------------------------------------------------------------------
//   writeState_
//-----------------------------------------------------------------------

// Write current (converged) step solution vector (STATE)
    

void ParaviewModule::writeState_ ( Ref<PrintWriter> vtuWriter,
                                           const Properties& globdat )
{   
   // Displacement field

   Vector  state;
   StateVector::get ( state, dofs_, globdat );

   int maxDofCount = rank_ == 2 ? rank_+1 : rank_;

   const String str1 =  String::format("<DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , maxDofCount);
   *vtuWriter          << str1;

   IdxVector  idofs    ( rank_       );
   Vector     nodeDisp ( maxDofCount );

   for (idx_t i = 0; i < nodeCount_; ++i)
   {
    nodeDisp = 0.;
    try
    {
      dofs_->getDofIndices ( idofs, i, dofTypes_[slice(BEGIN, rank_)] );
      nodeDisp[slice(BEGIN,rank_)] = select ( state, idofs );
    }
    catch ( const jem::IllegalInputException& ex )
    {}

    for (idx_t j = 0; j< maxDofCount; ++j)
    {
      *vtuWriter << String::format( "%12.6e   ", nodeDisp[j] );
    }
    *vtuWriter << endl;
   }

   *vtuWriter << "</DataArray> \n";


   // (Possible) additional fields
   // Field names chosen from FE Model

   if ( dofTypeCount_ > rank_ )
   {

    int addlDofs = dofTypeCount_ - rank_;

    for ( int j = 1; j <= addlDofs; j++ )
    {

      const String str2 = String::format("<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , dofNames_[rank_+j-1] , 1 );
      *vtuWriter << str2;

     idx_t    idof;
     double   val;

     for (idx_t i = 0; i < nodeCount_; ++i)
     {
       val = 0.0;
       try
       {
         idof = dofs_->getDofIndex ( i, dofTypes_[rank_+j-1] );
         val  = state[idof];
       }
       catch ( const jem::IllegalInputException& ex )
       {
         
       }
       *vtuWriter << String::format( "%12.6e   ", val );
       *vtuWriter << endl;
     }

     *vtuWriter << "</DataArray> \n";

    }
   }
}

//-----------------------------------------------------------------------
//   writeVTUClosure_
//-----------------------------------------------------------------------
    
void ParaviewModule::writeVTUClosure_ ( Ref<PrintWriter> vtuWriter )
{
   *vtuWriter << "</Piece> \n";
   *vtuWriter << "</UnstructuredGrid> \n";
   *vtuWriter << "</VTKFile> \n";
}

//-----------------------------------------------------------------------
//   writeTable_
//-----------------------------------------------------------------------
    

void ParaviewModule::writeTable_ 

     ( Ref<PrintWriter> vtuWriter,
       Ref<XTable>      table,
       const String&    tname )
{
   using jive::StringVector;

   StringVector     colNames = table->getColumnNames  ();
   idx_t            colCount = colNames.size          ();
   idx_t            rowCount = table->rowCount        ();

   Matrix           coord;

   // convert a table to a matrix (function in BasicUtils.cpp)
   
   matrixFromTable ( coord, *table, colNames );
   
   *vtuWriter << String::format(" <DataArray  type=\"Float64\"  Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n",
                             tname, colCount);

   for (idx_t i = 0; i < rowCount; ++i)
   {
      for (idx_t j = 0; j< colCount; ++j)
      *vtuWriter << String::format( "%12.8f   ", coord(i,j) );
      *vtuWriter << endl;
   }

   *vtuWriter << "</DataArray> \n";
}


//-----------------------------------------------------------------------
//   writeTable_ (overloaded)
//-----------------------------------------------------------------------
    

void ParaviewModule::writeTable_ 

     ( Ref<PrintWriter> vtuWriter,
       Ref<XTable>      table,
       const String&    tname,
       const int        rowCount )
{
   using jive::StringVector;

   StringVector     colNames = table->getColumnNames  ();
   idx_t            colCount = colNames.size          ();

   Matrix           coord;

   // convert a table to a matrix (function in BasicUtils.cpp)
   
   matrixFromTable ( coord, *table, colNames );
   
   *vtuWriter << String::format(" <DataArray  type=\"Float64\"  Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n",
                             tname, colCount);

   for (idx_t i = 0; i < rowCount; ++i)
   {
      for (idx_t j = 0; j< colCount; ++j)
      *vtuWriter << String::format( "%12.8f   ", coord(i,j) );
      *vtuWriter << endl;
   }

   *vtuWriter << "</DataArray> \n";
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void ParaviewModule::shutdown ( const Properties& globdat )
{
  // Do nothing!
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ParaviewModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  fileName_ = "problem_out";
  interval_ = 1;

  Properties  myProps = props.findProps ( myName_ );

  myProps.find ( fileName_,  FILENAME_PROP        );
  myProps.find ( interval_,  PRINT_INTERVAL_PROP  );
  myProps.find ( pointData_, POINT_DATA_PROP      );
  myProps.find ( cellData_,  CELL_DATA_PROP       );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ParaviewModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );
  
  myConf.set  ( FILENAME_PROP,        fileName_  );
  myConf.set  ( PRINT_INTERVAL_PROP,  interval_  );
  myConf.set  ( POINT_DATA_PROP,      pointData_ );
  myConf.set  ( CELL_DATA_PROP,       cellData_  );
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  ParaviewModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat )
{
  return newInstance<Self> ( name );
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareParaviewModule
//-----------------------------------------------------------------------

void declareParaviewModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( "paraview", & ParaviewModule::makeNew );
}
