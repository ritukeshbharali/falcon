
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

#include <utility>

/* Include user-defined class headers */

#include "ParaviewModule.h"
#include "util/TableUtils.h"
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

  const String context = getContext();

  // Get the MPI communicator
  mpx_    = Globdat::getMPContext ( globdat );

  // Declare the property tree for this module

  Properties  myProps = props.getProps ( myName_ );
  Properties  myConf  = conf.makeProps ( myName_ );

  // Get the element group and extract corresponding nodes and dof space

  egroup_ = ElemGroup::get      ( myConf, myProps, globdat, context );
  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();
  dofs_   = XDofSpace::get      ( nodes_.getData(), globdat );

  // Extract the relevant elements and nodes
  // inodes is a bit hacky (unavoidable)!

  ielems_  .resize( egroup_.size() );
  ielems_ = egroup_.getIndices();

  IdxVector inodes = elems_.getUniqueNodesOf 
                                   ( egroup_.getIndices() );
   
  inodes_.resize( inodes.size() );
  inodes_ = inodes;

  // Sort inodes in ascending order (without this,
  // the solution field is incorrectly mapped to nodes)

  jem::sort(inodes_);

  // Get some additional data ( number of nodes, elements
  // , mesh rank, dofs per node )
  
  nodeCount_    = inodes_.size     ();
  elemCount_    = ielems_.size     ();
  rank_         = nodes_.rank      ();
  dofTypeCount_ = dofs_->typeCount ();

  // Get nodes per element

  nodesPerElem_ = elems_.getElemNodeCount(ielems_[0]);

  // Get dof names and their types

   dofNames_.resize ( dofTypeCount_ );
   dofTypes_.resize ( dofTypeCount_ );
   
   dofNames_ = dofs_->getTypeNames();

   for ( idx_t i = 0; i < dofTypeCount_; i++ )
   {
     dofTypes_[i] = dofs_->findType ( dofNames_[i] );
   }

   // Element mapping based on mesh rank
   // (2D)

   if ( rank_ == 2 )
   {
     // Check for triangles (T3, T6)
     if ( nodesPerElem_ == 3 )
     {
       nVtxPerElem_  = 3;
       vtkCellCode_  = 5;

       vtxNodes_.resize( nVtxPerElem_ );
       vtxNodes_ = {0,1,2};
     }
     else if ( nodesPerElem_ == 6 )
     {
       nVtxPerElem_  = 3;
       vtkCellCode_  = 5;

       vtxNodes_.resize( nVtxPerElem_ );
       vtxNodes_ = {0,2,4};
     }

     // Check for Quads (Q4, Q8, Q9)
     else if ( nodesPerElem_ == 4 )
     {
       nVtxPerElem_  = 4;
       vtkCellCode_  = 9;

       vtxNodes_.resize( nVtxPerElem_ );
       vtxNodes_ = {0,1,2,3};
     }
     else if ( nodesPerElem_ == 8 || nodesPerElem_ == 9 )
     {
       nVtxPerElem_  = 4;
       vtkCellCode_  = 9;

       vtxNodes_.resize( nVtxPerElem_ );
       vtxNodes_ = {0,2,4,6};
     }
     else
     {
       const String  context = getContext ();

       throw IllegalInputException (
         context,
         String::format (
           "invalid 2D element (supported types: T3, T6, Q4, Q8, Q9)"
        )
       );
     }
   }
   // (3D)
   else if ( rank_ == 3 )
   {
    // Check for Tetrahedrons (Tet4, Tet10)
    if ( nodesPerElem_ == 4 )
    {
      nVtxPerElem_  = 4;
      vtkCellCode_  = 10;

      vtxNodes_.resize ( nVtxPerElem_ );
      vtxNodes_ = {0,1,2,3};
    }
    else if ( nodesPerElem_ == 10 )
    {
      nVtxPerElem_  = 4;
      vtkCellCode_  = 10;

      vtxNodes_.resize ( nVtxPerElem_ );
      vtxNodes_ = {0,2,4,9};
    }

    // Check for Hexahedrons (Hex8, Hex20)
    else if ( nodesPerElem_ == 8 )
    {
      nVtxPerElem_  = 8;
      vtkCellCode_  = 12;

      vtxNodes_.resize ( nVtxPerElem_ );
      vtxNodes_ = {0,1,2,3,4,5,6,7};
    }
    else if ( nodesPerElem_ == 20 )
    {
      nVtxPerElem_  = 8;
      vtkCellCode_  = 12;

      vtxNodes_.resize ( nVtxPerElem_ );
      vtxNodes_ = {0,2,4,6,12,14,16,18};
    }
    else
    {
      const String  context = getContext ();

      throw IllegalInputException (
        context,
        String::format (
           "invalid 3D element (supported types: Tet4, Tet10, Hex8, Hex20)"
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
  // Write a PVD file ( overwritten on every step )
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
  // Write a VTU file ( for current step )
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

  // Write PointData header. This is because, the solution
  // which is PointData (nodal data) is always printed to
  // the vtu file by default.

  *vtuWriter << "<PointData  Vectors=\"NodalData\"> \n";

  // Write solution (STATE0) to file. Instead of globdat,
  // send the state vector.

  writeState_ ( vtuWriter, globdat );

  // Proceed to write additional PointData requested
  // by the user.

  if ( pointData_.size() > 0 )
  {
    Properties      params  ("actionParams");
    Vector          weights ( nodeCount_ );

    // Loop over each point data requested

    for ( idx_t iData = 0; iData < pointData_.size(); iData++ )
    {
      String               tname ( pointData_[iData] );
      weights            = 0.0;

      Ref<Model>  model  = Model::get ( globdat, getContext() );
      Ref<XTable> xtable = newInstance<DenseTable>  ( tname, 
                                             nodes_.getData() );

      params.set ( ActionParams::TABLE,         xtable  );
      params.set ( ActionParams::TABLE_NAME,    tname   );
      params.set ( ActionParams::TABLE_WEIGHTS, weights );

      // Ask models to update the table

      bool updated = model->takeAction ( Actions::GET_TABLE, 
                                              params, globdat );

      // If table is updated, scale the rows and write to file

      if ( updated )
      {
        weights = where ( abs( weights ) < 
                               Limits<double>::TINY_VALUE,
                                           1.0, 1.0 / weights );
        
        // If dummy nodes are present, xtable will have more
        // rows than the size of weights vector. Consequently,
        // xtable -> scaleRows ( weights ) will return errors.
        // As a quick fix, we extend weights vector to match
        // the table size.

        if ( xtable->rowCount() > weights.size() )
        {
          weights.reshape( xtable->rowCount() );

          weights[slice(weights.size(),END)] = 0.0;
        }

        xtable -> scaleRows    ( weights );

        // Write table excluding dummy nodes

        writeTable_ ( vtuWriter, xtable, tname, nodeCount_ );
      }
    }  
  }

  // Write PointData closure

  *vtuWriter << "</PointData> \n";

  // Proceed to write additional CellData requested
  // by the user.

  if ( cellData_.size() > 0 )
  {
    Properties      params  ("actionParams");
    Vector          weights ( elemCount_   );

    // Write CellData header to file

    *vtuWriter << "<CellData> \n";

    // Loop over each cell data requested

    for ( idx_t iData = 0; iData < cellData_.size(); iData++ )
    {
      String               tname ( cellData_[iData] );
      weights            = 1.0;

      Ref<Model>  model  = Model::get ( globdat, getContext() );
      Ref<XTable> xtable = newInstance<DenseTable>  ( tname,
                                             elems_.getData() );

      params.set ( ActionParams::TABLE,         xtable  );
      params.set ( ActionParams::TABLE_NAME,    tname   );
      params.set ( ActionParams::TABLE_WEIGHTS, weights );

      bool updated = model->takeAction ( Actions::GET_TABLE,
                                              params, globdat );

      // If table is updated, scale the rows and write to file

      if ( updated )
      {
        //weights = where ( abs( weights ) < 
        //                       Limits<double>::TINY_VALUE,
        //                                   1.0, 1.0 / weights );
        //xtable -> scaleRows    ( weights );

        writeTable_   ( vtuWriter, xtable, tname, elemCount_ );
      }
    }

    // Write CellData closure to file

    *vtuWriter << "</CellData> \n";
  }

  // TO-DO: Add support for tensor data 
  // (e.g. for stress, strain)

  // Write closure for the VTU file

  writeVTUClosure_   ( vtuWriter );

  // Flush and Close the VTU file

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
  // Get nodal coordinates (includes dummy nodes too!)
  Matrix coords ( rank_, nodes_.size() );
  nodes_.getSomeCoords ( coords, inodes_ );

  // Write nodal data

  *vtuWriter << "<Points> \n";
  *vtuWriter << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\" format=\"ascii\" >\n";

  // Mesh rank dependent print
  // (2D)
  if ( rank_ == 2 )
  {
    for ( idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuWriter << coords(0,i) << " " 
                 << coords(1,i) << " " 
                 << 0. << '\n';
    }
  }
  // (3D)
  else if ( rank_ == 3 )
  {
    for ( idx_t i = 0; i < nodeCount_; ++i )
    {
      *vtuWriter << coords(0,i) << " " 
                 << coords(1,i) << " " 
                 << coords(2,i) << '\n';
    }
  }

  *vtuWriter << "</DataArray>\n </Points>\n";

  // Initialize inodes to store element connectivity

  IdxVector inodes(nodesPerElem_);

  // Write element connectivity header

  *vtuWriter << "<Cells> \n";
  *vtuWriter << "<DataArray  type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

  // Write element connectivity data
  // (Only the vertices)
  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    const int ielem = ielems_[ie];

    elems_.getElemNodes ( inodes, ielem );

    for ( idx_t i = 0; i < nVtxPerElem_; i++ )
    {
      *vtuWriter << inodes[vtxNodes_[i]] << " ";
    }

    *vtuWriter << endl;
  }

  *vtuWriter << "</DataArray>\n";

  // Write cell (element) offset

  idx_t offset = 0;

  *vtuWriter << "<DataArray  type=\"Int32\"  Name=\"offsets\"  format=\"ascii\"> \n";

  for ( idx_t ie = 0; ie < elemCount_; ie++ )
  {
    offset += nVtxPerElem_;
    *vtuWriter << offset << endl;
  }

  *vtuWriter << "</DataArray>\n";

  // Write cell (element) types [VTK notation]

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

// Writes the current (converged) step solution vector (STATE)

void ParaviewModule::writeState_ ( Ref<PrintWriter> vtuWriter,
                              const Properties& globdat  )
{
  // Get the state vector

  Vector state; StateVector::get ( state, dofs_, globdat );

  // Allocate vectors to extract dof number and values

  IdxVector idof;  idof .resize( 1 );
  IdxVector idofs; idofs.resize( nodeCount_ );
  Vector    ivals; ivals.resize( nodeCount_ );

  // Loop over the dofs

  for ( idx_t i = 0; i < dofTypeCount_; i++ )
  {
    // Dof type index

    idof[0] = i;

    // Get the dofs

    try
    {
      dofs_->getDofIndices ( idofs, inodes_, idof );
      ivals = select ( state, idofs );

      // Write header to file

      const String headStr = String::format("<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\"> \n" , dofNames_[i] , 1 );
      *vtuWriter << headStr;

      // Write the data

      for (idx_t j = 0; j< nodeCount_-1; j++)
      {
        *vtuWriter << String::format( "%12.6e  \n", ivals[j] );
      }
      *vtuWriter << String::format( "%12.6e  ", ivals[nodeCount_-1] );
      *vtuWriter << endl;

      // Write closure

      *vtuWriter << "</DataArray> \n";

    }
    catch ( const jem::IllegalInputException& ex )
    {}
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
      {
        *vtuWriter << String::format( "%12.8f   ", coord(i,j) );
      }
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

  ModuleFactory::declare ( "Paraview", & ParaviewModule::makeNew );
}
