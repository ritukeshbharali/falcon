
/** @file GradientCrystalPlasticityROM.cpp
 *  @brief Crystal visco-plasticity reduced order model with gradient regularization.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 25 July 2024
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2024] (RB)
 * 
 *  @todo Add post-processing options (displacement, strains, stress)
 *
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/Float.h>
#include <jem/base/IllegalInputException.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/util/Timer.h>
#include <jive/model/Actions.h>

/* Include other headers */

#include "FalconSolidMechROMs.h"
#include "GradientCrystalPlasticityROM.h"


//=======================================================================
//   class GradientCrystalPlasticityROM
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  GradientCrystalPlasticityROM::DISP_NAMES[3]   = { "dx", "dy", "dz" };

const char*  GradientCrystalPlasticityROM::SHAPE_PROP      = "shape";
const char*  GradientCrystalPlasticityROM::DTIME_PROP      = "dtime";
const char*  GradientCrystalPlasticityROM::NSLIPS_PROP     = "nslips";
const char*  GradientCrystalPlasticityROM::NMODES_PROP     = "nmodes";
const char*  GradientCrystalPlasticityROM::MODE_NODES_PROP = "modeNodes";
const char*  GradientCrystalPlasticityROM::SLIP_HAT0_FILE  = "slipHat0File";
const char*  GradientCrystalPlasticityROM::BASIS           = "basis";
const char*  GradientCrystalPlasticityROM::SLIP_HAT_FILE   = "slipHatFile";
const char*  GradientCrystalPlasticityROM::TAU_HAT_FILE    = "tauHatFile";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


GradientCrystalPlasticityROM::GradientCrystalPlasticityROM

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super      ( name ),
    nslips_    ( 1 ),
    dtime_     ( 1.0 ),
    nmodes_    ( 1 )

{
  using jem::util::Timer;
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );

  const String  context = getContext ();

  // Get the element group assigned to this model.

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();
  rank_   = nodes_ .rank        ();

  // Make sure that the number of spatial dimensions (the rank of the
  // mesh) is valid.

  if ( rank_ < 1 || rank_ > 3 )
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  }

  // Create an internal shape object for computing the element shape
  // functions.

  shape_ = IShapeFactory::newInstance (
    joinNames ( myName_, SHAPE_PROP ),
    conf,
    props
  );

  // Make sure that the rank of the shape matches the rank of the
  // mesh.

  if ( shape_->globalRank() != rank_ )
  {
    throw IllegalInputException (
      context,
      String::format (
        "shape has invalid rank: %d (should be %d)",
        shape_->globalRank (),
        rank_
      )
    );
  }

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements (
    context,
    egroup_.getIndices (),
    shape_->nodeCount  ()
  );

  // Compute the total number of integration points.

  const int ipCount = shape_->ipointCount() * egroup_.size();

  // Get the number of slip systems

  myProps.get ( nslips_, NSLIPS_PROP );
  myConf .set ( NSLIPS_PROP, nslips_ );

  // Get the time step-size

  myProps.get ( dtime_,  DTIME_PROP );
  myConf. set ( DTIME_PROP, dtime_  );

  // Get some material properties

  myProps.get ( tstar_, "tstar" );
  myProps.get ( tauY_,  "tauY"  );
  myProps.get ( n_,     "n"     );

  myConf.set ( "tstar", tstar_ );
  myConf.set ( "tauY",  tauY_  );
  myConf.set ( "n",     n_     );

  // Get the number of POD modes 

  myProps.get ( nmodes_, NMODES_PROP );
  myConf. set ( NMODES_PROP, nmodes_ );

  // Get the dummy node group for modes

  myProps.get ( modeNGroup_, MODE_NODES_PROP );

  modeNodes_  = NodeGroup::get ( modeNGroup_, nodes_, 
                                   globdat, context );

  // Check whether node group size matches number of modes

  if ( modeNodes_.size() != nmodes_ )
  {
    throw IllegalInputException (
      context,
      "wrong modeNodes size!"
    );
  }

  // Get the XDofSpace

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  // dofTypes has one element, the mode amplitude

  dofTypes_.resize( 1 );

  // Make a shallow copy of dofTypes_ for mode dofs
  // (Maybe not required!)

  modeTypes_.ref ( dofTypes_[slice(BEGIN, 1)] );

  // Assign a dof name and add dof type to XDofSpace

  String dofName = "mode";
  modeTypes_[0]  = dofs_->addType ( dofName );

  // Assign dofs to the dummy mode nodes

  dofs_->addDofs ( modeNodes_.getIndices(),
                   dofTypes_[slice(BEGIN, 1)] 
                 );

  myProps.get ( slipHat0File_, SLIP_HAT0_FILE );

  // Now, the pre-computed basis from offline data
  // are parsed and loaded onto Jive Matrix, Cubix
  // data structures.

  Timer      t;  t.start();

  // Parse gamma hat 0 (basis corresponding to unit
  // load)

  slipHat0_.resize( ipCount, nslips_ );

  System::out() << "Parsing " << slipHat0File_ << "... \n";

  // Open file, return error if it fails
  iFile_ = newInstance<FileReader> ( slipHat0File_ );

  // Parse entries
  for ( int irow = 0; irow < ipCount; irow++ )
  {
    for ( int icol = 0; icol < nslips_; icol++ )
    {
      slipHat0_(irow,icol) = iFile_->parseFloat();
    }
  }

  // Close file
  iFile_->close();

  myConf.set ( SLIP_HAT0_FILE, slipHat0File_ );

  // Parse all gamma hat sensitivities from file

  myProps.get ( slipHatFile_,  SLIP_HAT_FILE  );
  myProps.get ( basis_, BASIS );

  // Check whether node group size matches number of modes

  if ( basis_.size() != nmodes_ )
  {
    throw IllegalInputException (
      context,
      "size of basis does not match number of modes!"
    );
  }

  myConf.set ( BASIS, basis_ );

  slipHat_.resize( ipCount, nslips_, nmodes_ );

  for ( int imode = 0; imode < nmodes_; imode++ )
  {
    // Add mode index to file name
    String fileName = slipHatFile_ + String(basis_[imode]);

    System::out() << "Parsing " << fileName << "... \n";

    // Open file, return error if it fails
    iFile_ = newInstance<FileReader> ( fileName );

    // Parse entries
    for ( int irow = 0; irow < ipCount; irow++ )
    {
      for ( int icol = 0; icol < nslips_; icol++ )
      {
        slipHat_(irow,icol,imode) = iFile_->parseFloat();
      }
    }

    // Close file
    iFile_->close();
  }

  myConf.set ( SLIP_HAT_FILE,  slipHatFile_ );

  // Parse all tau hat sensitivities from file

  myProps.get ( tauHatFile_,  TAU_HAT_FILE  );

  tauHat_.resize( ipCount, nslips_, nmodes_ );

  for ( int imode = 0; imode < nmodes_; imode++ )
  {
    // Add mode index to file name
    String fileName = tauHatFile_ + String(basis_[imode]);

    System::out() << "Parsing " << fileName << "... \n";

    // Open file, return error if it fails
    iFile_ = newInstance<FileReader> ( fileName );

    // Parse entries
    for ( int irow = 0; irow < ipCount; irow++ )
    {
      for ( int icol = 0; icol < nslips_; icol++ )
      {
        tauHat_(irow,icol,imode) = iFile_->parseFloat();
      }
    }

    // Close file
    iFile_->close();
  }

  myConf.set ( TAU_HAT_FILE, tauHatFile_ );

  System::out() << "Parsed pre-computed data in " << t.toDouble()
                << " seconds ...\n";
}


GradientCrystalPlasticityROM::~GradientCrystalPlasticityROM ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void GradientCrystalPlasticityROM::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void GradientCrystalPlasticityROM::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool GradientCrystalPlasticityROM::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // Compute and store fixed parts of jacobian matrix
  // and residual during init stage

  if ( action == Actions::INIT )
  {
    init_ ();

    return true;
  }

  // Compute the internal force vector

  if ( action == Actions::GET_INT_VECTOR )
  {
    Vector  state, state0;
    Vector  force; 

    // Get the current and previous step state.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );
    
    // Get the internal force vector.

    params.get ( force, ActionParams::INT_VECTOR );

    getMatrix_ ( nullptr, force, state, state0 );

    return true;
  }

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {
    Ref<MatrixBuilder>  mbuilder;

    Vector  state, state0;
    Vector  force; 

    // Get the current and previous step state.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );
    
    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( mbuilder, force, state, state0 );

    return true;
  }

  // compute mass matrix 

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );
    
    getMatrix2_( *mbuilder );

    return true;
  }

  /** COMMIT: Requests an action when a computation step has converged.
   *  Ideally, at this point, the material internal state variables
   *  are swapped.
  */ 

  if ( action == Actions::COMMIT )
  {
    Vector  state;

    // Get the current state.

    StateVector::get ( state,  dofs_, globdat );

    System::out() << "solution: " << state << "\n";

    step_++;

    // Error

    errTotal_ += dtime_ * err_;

    System::out() << "error: " << errTotal_ << "\n";

    return true;
  }

  /** CHECK_COMMIT: Can used to discard the current step
  */ 

  if ( action == Actions::CHECK_COMMIT )
  {
    return true;
  }

  /** GET_TABLE: Requests a post-processing action, such as computing
   *  tabular data (stress, strains, etc.) to print VTK files.
  */ 
  
  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  return false;
}


//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void GradientCrystalPlasticityROM::init_ ()

{
  // Re-size matrix M_ and vector fext_
  M_   .resize( nmodes_, nmodes_ ); M_    = 0.0;
  fext_.resize( nmodes_ );          fext_ = 0.0;           

  // Get the elements associated with this model
  IdxVector   ielems     = egroup_ .getIndices ();

  // Get the mode nodes for this model
  IdxVector   modenodes  = modeNodes_.getIndices ();

  // Get some additional details
  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   mNodeCount = modenodes.size      ();

  // Extract the solution (probably not required!)

  IdxVector   modeDofs ( mNodeCount );
  dofs_->getDofIndices ( modeDofs, modenodes, modeTypes_ );

  // Allocate inodes, coords, ipWeights, wip
  IdxVector   inodes     ( nodeCount );
  Matrix      coords     ( rank_, nodeCount );
  Vector      ipWeights  ( ipCount   );

  // Iterate over all elements assigned to this model

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    // Get the element coordinates.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );

    // Get the element ip weights
    shape_-> getIntegrationWeights( ipWeights, coords );

    // Compute M_
    for ( int islip = 0; islip < nslips_; islip++ )
    {
      M_    += matmul( tauHat_  (ielem,islip,ALL), 
                    slipHat_(ielem,islip,ALL) )
               * ipWeights[0];

      fext_ += tauHat_  (ielem,islip,ALL) 
               * slipHat0_(ielem,islip)
               * ipWeights[0];
    }
  }

  // Scale with dtime
  M_ /= dtime_;

  // Initialize step counter and error

  step_     = 0;
  err_      = 0.0;
  errTotal_ = 0.0;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void GradientCrystalPlasticityROM::getMatrix_

  ( Ref<MatrixBuilder>  mbuilder,
    const Vector&       force,
    const Vector&       state,
    const Vector&       state0 )

{
  // Get the elements associated with this model
  IdxVector   ielems     = egroup_ .getIndices ();

  // Get the mode nodes for this model
  IdxVector   modenodes  = modeNodes_.getIndices ();

  // Get some additional details
  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   mNodeCount = modenodes.size      ();

  IdxVector   modeDofs ( mNodeCount );
  dofs_->getDofIndices ( modeDofs, modenodes, modeTypes_ );

  // Allocate inodes, coords, ipWeights, wip
  IdxVector   inodes     ( nodeCount );
  Matrix      coords     ( rank_, nodeCount );
  Vector      ipWeights  ( ipCount   );

  MChain1  mc1;

  Matrix   mat ( mNodeCount, mNodeCount ); 

  mat   = M_;
  force = 0.0;
  err_  = 0.0;

  // Iterate over all elements assigned to this model

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    // Get the element coordinates.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );

    // Get the element ip weights
    shape_-> getIntegrationWeights( ipWeights, coords );

    // Allocate tau, phi and dphi

    Vector tau  ( nslips_ ); tau  = 0.0;
    Vector phi  ( nslips_ ); phi  = 0.0;
    Vector dphi ( nslips_ ); dphi = 0.0;

    // Compute tau, phi and dphi

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      tau[islip]  = dot( tauHat_(ielem, islip,ALL), state );
      
      phi[islip]  = sign( tau[islip] ) / tstar_
                    * ::pow( fabs( tau[islip] ) / 
                              tauY_, n_ );

      dphi[islip] = 1.0 / tstar_ * n_ / tauY_
                    * ::pow( fabs( tau[islip] ) / 
                             tauY_, n_ - 1.0 );

      force      -= ipWeights[0] * tauHat_(ielem,islip,ALL)
                    * phi[islip];

      mat        -= matmul( tauHat_ (ielem,islip,ALL), 
                    slipHat_(ielem,islip,ALL) )
                    * ipWeights[0] * dphi[islip];

      // Compute slip and slip0 (required for error)

      double slip  = slipHat0_(ielem,islip) * (step_+1) * dtime_
                    + dot( slipHat_(ielem, islip,ALL), state );
      double slip0 = slipHat0_(ielem,islip) * step_ * dtime_
                    + dot( slipHat_(ielem, islip,ALL), state0 );

      // Error computation

      err_ += ipWeights[0] * 
             ::pow( (( slip - slip0 )/dtime_ - phi[islip]), 2.0 )
             / (dphi[islip] + 1.e-06);
    }
  }

  // Update force with pre-computed parts

  force += mc1.matmul( M_, Vector(state - state0) )
           +  fext_;

  // Assemble global

  if ( mbuilder != nullptr )
  {
    mbuilder -> addBlock ( modeDofs, modeDofs, mat );
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void GradientCrystalPlasticityROM::getMatrix2_

    ( MatrixBuilder&          mbuilder )
{
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool GradientCrystalPlasticityROM::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  Ref<XTable>  table;
  Vector       weights;
  String       name;


  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table,   ActionParams::TABLE );
  params.get ( name,    ActionParams::TABLE_NAME );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );

  // Element tau
  if ( name == "tau" && 
       table->getRowItems() == elems_.getData() )
  {
    Vector  state;
    StateVector::get ( state, dofs_, globdat  );

    getElemTau_ ( *table, weights, state );

    return true;
  }

  // Element slip
  if ( name == "slip" && 
       table->getRowItems() == elems_.getData() )
  {
    Vector  state;
    StateVector::get ( state, dofs_, globdat  );

    getElemSlip_ ( *table, weights, state );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getElemTau_
//-----------------------------------------------------------------------


void GradientCrystalPlasticityROM::getElemTau_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  state )

{
  // Get the elements associated with this model
  IdxVector  ielems     = egroup_.getIndices  ();

  // Get the mode nodes for this model
  IdxVector   modenodes  = modeNodes_.getIndices ();

  // Get some additional details
  const int   elemCount  = ielems.size         ();
  //const int   mNodeCount = modenodes.size      ();

  // Extract the solution (probably not required!)

  //IdxVector   modeDofs ( mNodeCount );
  //dofs_->getDofIndices ( modeDofs, modenodes, modeTypes_ );

  //Vector sol  ( mNodeCount ); sol  = select ( state,  modeDofs );

  Matrix     elTau      ( elemCount, nslips_ );
  IdxVector  jcols      ( nslips_   );

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    String tauName = "tau" + String(islip);
    jcols[islip]   = table.addColumn ( tauName );
  }

  MChain1    mc1;

  elTau = 0.0;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      // Assuming one point Gauss integration

      elTau(ielem,islip) += dot( tauHat_(ielem, islip,ALL), state );
    }
  }

  // Add the tau to the table.

  table.setBlock ( ielems, jcols, elTau );
}


//-----------------------------------------------------------------------
//   getElemSlip_
//-----------------------------------------------------------------------


void GradientCrystalPlasticityROM::getElemSlip_

  ( XTable&        table,
    const Vector&  weights,
    const Vector&  state )

{
  // Get the elements associated with this model
  IdxVector  ielems     = egroup_.getIndices  ();

  // Get the mode nodes for this model
  IdxVector   modenodes  = modeNodes_.getIndices ();

  // Get some additional details
  const int   elemCount  = ielems.size         ();

  // Initialize 

  Matrix     elSlip     ( elemCount, nslips_ );
  IdxVector  jcols      ( nslips_   );

  for ( int islip = 0; islip < nslips_; islip++ )
  {
    String slipName = "slip" + String(islip);
    jcols[islip]    = table.addColumn ( slipName );
  }

  MChain1    mc1;

  elSlip = 0.0;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < elemCount; ie++ )
  {
    // Get the global element index.

    const int  ielem = ielems[ie];

    for ( int islip = 0; islip < nslips_; islip++ )
    {
      // Assuming one point Gauss integration
      // A scaling dtime_ is included for slipHat0

      elSlip(ielem,islip) += slipHat0_(ielem,islip) * step_ * dtime_
                            + dot( slipHat_(ielem, islip,ALL), state );
    }
  }

  // Add the slips to the table.

  table.setBlock ( ielems, jcols, elSlip );
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void GradientCrystalPlasticityROM::checkCommit_

  ( const Properties&  params )

{
  // System::info() << myName_ << " : check commit ... do nothing!\n";
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newGradientCrystalPlasticityROM
//-----------------------------------------------------------------------


static Ref<Model>     newGradientCrystalPlasticityROM

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<GradientCrystalPlasticityROM> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareGradientCrystalPlasticityROM
//-----------------------------------------------------------------------


void declareGradientCrystalPlasticityROM ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "GradientCrystalPlasticityROM", & newGradientCrystalPlasticityROM );
}
