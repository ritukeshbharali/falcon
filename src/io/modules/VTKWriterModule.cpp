
/** @file VTKWriterModule.cpp
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

#include "VTKWriterModule.h"
#include "util/BasicUtils.h"

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
//   class VTKWriterModule
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  VTKWriterModule::FILE_PROP         = "fileName";
const char*  VTKWriterModule::INTERVAL_PROP     = "interval";
const char*  VTKWriterModule::DATA_PROP         = "data";
const char*  VTKWriterModule::DATA_TYPE_PROP    = "dataType";



//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


VTKWriterModule::VTKWriterModule ( const String& name ) :
  Super ( name )
{}

VTKWriterModule::~VTKWriterModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status VTKWriterModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
   using jive::mp::Globdat;

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

   const StringVector dofNames = dofs_->getTypeNames();

   dofNames_.resize ( dofTypeCount_ );
   dofTypes_.resize ( dofTypeCount_ );
  
   dofNames_ = dofNames;

   for ( idx_t i = 0; i < dofTypeCount_; i++ )
   {
     dofTypes_[i] = dofs_->addType ( dofNames_[i] );
   }

   // MPI ( Get the number of process )  

   mpx_    = Globdat::getMPContext ( globdat );
   nProcs_ = 1;
   nProcs_ = mpx_-> size ();

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
       cnodes_.resize( numVertexesPerCell_ );
       cnodes_ = {0,2,4};
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
       cnodes_.resize( numVertexesPerCell_ );
       cnodes_ = {0,2,4,6};
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
       cnodes_.resize( numVertexesPerCell_ );
       cnodes_ = {0,2,4,9}; 
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
       cnodes_.resize( numVertexesPerCell_ );
       cnodes_ = {0,2,4,6,12,14,16,18};
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

   // Vector to store all timestamps

   timeStamps_.resize(1);

   return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------

Module::Status VTKWriterModule::run ( const Properties& globdat )
{
  using jem::Ref;
  using jem::System;
  using jive::Vector;
  using jive::util::ItemSet;
  using jive::model::Model;
  using jive::model::ActionParams;
  using jive::model::Actions;

  const String   context = getContext ();

  // Get current process

  int proc                   = mpx_-> myRank ();

  // Get current (converged) step and time

  int           step;
  double        time;

           globdat.get  ( step, Globdat::TIME_STEP );
  bool t = globdat.find ( time,  Globdat::TIME     );

  if ( t )
  {
    timeStamps_.pushBack( time );
  }
  else
  {
    timeStamps_.pushBack( step );
  }

  // Print first step and then every 'interval_' steps

  if ( ( step % interval_ ) != 0 && step > 1 ) return OK;

  const String fileName = fileName_ + "_p" + String(proc) 
                               + "_" + String(step)  + ".vtu";

  // Start a PVD writer (for timestamps)

  Ref<PrintWriter>         pvdFile_;
  pvdFile_  = newInstance<PrintWriter> (
              newInstance<FileWriter>  ( 
                           fileName_ + ".pvd" ));
  *pvdFile_ << "<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n";
  *pvdFile_ << "<Collection>\n";

  for (idx_t i = 0; i <= step; ++i)
   {
      if ( ( i % interval_ ) != 0 && step > 1 ) continue;  

      String fileNamei = fileName_ + "_p" + String(proc) 
                               + "_" + String(i)  + ".vtu";                             

      *pvdFile_ << "<DataSet file='"+ fileNamei 
               + "' groups='' part='0' timestep='"
               + String(timeStamps_[i])+"'/>\n";
    }

    // Close the pvd file
    *pvdFile_ << "</Collection>\n";
    *pvdFile_ << "</VTKFile>\n";

    pvdFile_->flush();
    pvdFile_->close();

  // Create a vtu writer

  Ref<PrintWriter> vtuFile   = newInstance<PrintWriter> ( 
                               newInstance<FileWriter>  ( 
                               fileName ) );

  // Write header, mesh and solution (STATE)

  writeHeader_ ( vtuFile );
  writeMesh_   ( vtuFile );
  writeState_  ( vtuFile, globdat );

  //System::out() << "Done printing state!! \n";

  // Write PointData (Nodal average for Gauss point variables)
  if ( dataType_ == "nodes" )
  {

    /** <PointData  Vectors="NodalData"> not required since
     *  it has been written by writeState_. However,
     *  PointData must be closed with </PointData> before
     *  moving on to <CellData>
    */

    const idx_t dataCount = dataNames_.size();

    Properties      params  ("actionParams");
    Vector          weights ( nodeCount_ );

      for ( idx_t i = 0; i < dataCount; i++ )
      {

        String              tname (dataNames_[i]);
        weights           = 0.0;
        Ref<Model> model  = Model::get ( globdat, getContext() );

        Ref<XTable>     xtable = newInstance<DenseTable>  ( tname, nodes_.getData() );

        params.set ( ActionParams::TABLE,         xtable  );
        params.set ( ActionParams::TABLE_NAME,    tname   );
        params.set ( ActionParams::TABLE_WEIGHTS, weights );

        bool updated = model->takeAction ( Actions::GET_TABLE, params, globdat );
      
        weights = where ( abs( weights ) < Limits<double>::TINY_VALUE, 1.0, 1.0 / weights );
        xtable -> scaleRows    ( weights );

        if ( updated ) 
        {
          writeTable_   ( vtuFile, xtable, tname );
        }
      }

    // Close PointData
    *vtuFile << "</PointData> \n";

  }

  // Write CellData (Gauss point average)
  else if ( dataType_ == "elems" )
  {

    // Close PointData started by writeState_
    *vtuFile << "</PointData> \n";

    // Start Cell Data
    *vtuFile << "<CellData> \n";

    const idx_t dataCount = dataNames_.size();

    Properties      params  ("actionParams");
    Vector          weights ( elemCount_ );

    for ( idx_t i = 0; i < dataCount; i++ )
      {

        String              tname (dataNames_[i]);
        weights           = 1.0;
        Ref<Model> model  = Model::get ( globdat, getContext() );

        Ref<XTable>     xtable = newInstance<DenseTable>  ( tname, elems_.getData() );

        params.set ( ActionParams::TABLE,         xtable  );
        params.set ( ActionParams::TABLE_NAME,    tname   );
        params.set ( ActionParams::TABLE_WEIGHTS, weights );

        bool updated = model->takeAction ( Actions::GET_TABLE, params, globdat );

        if ( updated ) 
        {
          writeTable_   ( vtuFile, xtable, tname, elemCount_ );
        }
      }

    // Close CellData
    *vtuFile << "</CellData> \n";
  }

  // Write Closure

  writeClosure_   ( vtuFile );

  // Flush and Close file

  vtuFile->flush();
  vtuFile->close();

  return OK;

}


//-----------------------------------------------------------------------
//   writeHeader_
//-----------------------------------------------------------------------
    

void VTKWriterModule::writeHeader_ ( Ref<PrintWriter> vtuFile )
{
  const String str1 = String::format("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"> \n"
                                 "<UnstructuredGrid> \n"
                                   "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\"> \n", 
                             nodeCount_, elemCount_ );

  *vtuFile << str1;
}



//-----------------------------------------------------------------------
//   writeMesh_
//-----------------------------------------------------------------------
    

void VTKWriterModule::writeMesh_ ( Ref<PrintWriter> vtuFile )
{

  // Write node coordinates

  *vtuFile << "<Points> \n";
  *vtuFile << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" >\n";

  Matrix coords(rank_,nodeCount_);

  nodes_.getCoords ( coords );

  if ( rank_ == 2 )
  {
    for(idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuFile << coords(0,i) << " " << coords(1,i) << " " << 0. << '\n';
    }
  }
  else
  {
    for(idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuFile << coords(0,i) << " " << coords(1,i) << " " << coords(2,i) << '\n';
    }
  }

  *vtuFile << "</DataArray>\n </Points>\n";

  // Write element connectivity

  IdxVector inodes(noNodePerElem_);

  *vtuFile << "<Cells> \n";
  *vtuFile << "<DataArray  type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  { 
    int  ielem = ielems_[ie]; 
    elems_.getElemNodes  ( inodes, ielem );

    // Write only the corner nodes for higher order elements
    // Middle nodes are not required for mesh visualization

    if (noNodePerElem_ == numVertexesPerCell_ )
    {
      for ( idx_t in = 0; in < noNodePerElem_; in++ ){
        *vtuFile << inodes[in] << " "; 
      }
    }
    else
    {
      for ( idx_t in = 0; in < numVertexesPerCell_; in++ ){
        *vtuFile << inodes[cnodes_[in]] << " "; 
      }
    }

    *vtuFile << endl;
  }

  *vtuFile << "</DataArray>\n";

  // write cell offset

  *vtuFile << "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\"> \n";

  int offset = 0;
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    offset += numVertexesPerCell_;
    *vtuFile <<  offset << endl;
  }

  *vtuFile << "</DataArray>\n";

  // Print cell types
  
  *vtuFile << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\"> \n";
  
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
     *vtuFile <<  vtkCellCode_ << endl;
  }

  *vtuFile << "</DataArray> \n </Cells> \n";

}

//-----------------------------------------------------------------------
//   writeState_
//-----------------------------------------------------------------------

// Write current (converged) step solution vector (STATE)
    

void VTKWriterModule::writeState_ ( Ref<PrintWriter> vtuFile,
                                           const Properties& globdat )
{
   *vtuFile << "<PointData  Vectors=\"NodalData\"> \n";
   
   // Displacement field

   Vector  state;
   StateVector::get ( state, dofs_, globdat );

   int maxDofCount = rank_ == 2 ? rank_+1 : rank_;

   const String str1 =  String::format("<DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , maxDofCount);
   *vtuFile          << str1;

   IdxVector  idofs    ( rank_       );
   Vector     nodeDisp ( maxDofCount );

   for (idx_t i = 0; i < nodeCount_; ++i)
   {
      nodeDisp = 0.;
      try{
        dofs_->getDofIndices ( idofs, i, dofTypes_[slice(BEGIN, rank_)] );

        nodeDisp[slice(BEGIN,rank_)] = select ( state, idofs );
      }
      catch ( const jem::IllegalInputException& ex )
      {
      }

      for (idx_t j = 0; j< maxDofCount; ++j)
      {
         *vtuFile << String::format( "%12.6e   ", nodeDisp[j] );
      }
      *vtuFile << endl;
   }

   *vtuFile << "</DataArray> \n";


   // (Possible) additional fields
   // Field names chosen from FE Model

   if ( dofTypeCount_ > rank_ )
   {

    int addlDofs = dofTypeCount_ - rank_;

    for ( int j = 1; j <= addlDofs; j++ )
    {

      const String str2 = String::format("<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , dofNames_[rank_+j-1] , 1 );
      *vtuFile << str2;

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
       *vtuFile << String::format( "%12.6e   ", val );
       *vtuFile << endl;
     }

     *vtuFile << "</DataArray> \n";

    }
   }
}

//-----------------------------------------------------------------------
//   writeClosure_
//-----------------------------------------------------------------------
    
void VTKWriterModule::writeClosure_ ( Ref<PrintWriter> vtuFile )
{
   *vtuFile << "</Piece> \n";
   *vtuFile << "</UnstructuredGrid> \n";
   *vtuFile << "</VTKFile> \n";
}

//-----------------------------------------------------------------------
//   writeTable_
//-----------------------------------------------------------------------
    

void VTKWriterModule::writeTable_ 

     ( Ref<PrintWriter> vtuFile,
       Ref<XTable>      table,
       const String&    tname )
{
   using jive::StringVector;

   StringVector     colNames = table->getColumnNames  ();
   idx_t            colCount = colNames.size          ();
   idx_t            rowCount = table->rowCount        ();

   Matrix           coord;

   // convert a table to a matrix (function in utilities.cpp)
   
   matrixFromTable ( coord, *table, colNames );
   
   *vtuFile << String::format(" <DataArray  type=\"Float64\"  Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n",
                             tname, colCount);

   for (idx_t i = 0; i < rowCount; ++i)
   {
      for (idx_t j = 0; j< colCount; ++j)
      *vtuFile << String::format( "%12.8f   ", coord(i,j) );
      *vtuFile << endl;
   }

   *vtuFile << "</DataArray> \n";
}


//-----------------------------------------------------------------------
//   writeTable_ (overloaded)
//-----------------------------------------------------------------------
    

void VTKWriterModule::writeTable_ 

     ( Ref<PrintWriter> vtuFile,
       Ref<XTable>      table,
       const String&    tname,
       const int        rowCount )
{
   using jive::StringVector;

   StringVector     colNames = table->getColumnNames  ();
   idx_t            colCount = colNames.size          ();

   Matrix           coord;

   // convert a table to a matrix (function in utilities.cpp)
   
   matrixFromTable ( coord, *table, colNames );
   
   *vtuFile << String::format(" <DataArray  type=\"Float64\"  Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n",
                             tname, colCount);

   for (idx_t i = 0; i < rowCount; ++i)
   {
      for (idx_t j = 0; j< colCount; ++j)
      *vtuFile << String::format( "%12.8f   ", coord(i,j) );
      *vtuFile << endl;
   }

   *vtuFile << "</DataArray> \n";
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void VTKWriterModule::shutdown ( const Properties& globdat )
{
  // Do nothing!
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void VTKWriterModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  interval_ = 1;
  dataType_ = "nodes";

  Properties  myProps = props.findProps ( myName_ );

  myProps.get  ( fileName_,  FILE_PROP      );
  myProps.find ( interval_,  INTERVAL_PROP  );
  myProps.get  ( dataNames_, DATA_PROP      );
  myProps.find ( dataType_,  DATA_TYPE_PROP );

  if ( dataType_ != "nodes" && dataType_ != "elems" )
  {
    const String  context = getContext ();

    throw IllegalInputException (
      context,
      String::format (
        "invalid data type: %s (should be nodes or elems)", dataType_
      )
    );

  }

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void VTKWriterModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );
  
  myConf.set  ( FILE_PROP,      fileName_  );
  myConf.set  ( INTERVAL_PROP,  interval_  );
  myConf.set  ( DATA_PROP,      dataNames_ );
  myConf.set  ( DATA_TYPE_PROP, dataType_  );
}

Ref<Module>  VTKWriterModule::makeNew

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
//   declareVTKWriterModule
//-----------------------------------------------------------------------

void declareVTKWriterModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( "vtkWriter", & VTKWriterModule::makeNew );
}
