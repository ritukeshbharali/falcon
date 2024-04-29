
/** @file AMGXSolver.cpp
 *  @brief Wrapper to solve a linear system with AMGX.
 * 
 *  Requires: 
 *     - NVIDIA GPU with Compute Capability >=3.0
 *     - AMGX installation with dependencies
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 02 January 2024 
 *
 *  Updates (when, what and who)
 *     - [XX YYYYY 2024]
 *       
 */

#if defined(WITH_AMGX)

/* Include jem and jive headers */

#include <cmath>
#include <jem/base/Array.h>
#include <jem/base/Float.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/ArithmeticException.h>
#include <jem/base/IllegalArgumentException.h>
#include <jem/base/IllegalOperationException.h>
#include <jive/util/error.h>
#include <jive/util/utilities.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/SparseMatrixObject.h>
#include <jive/solver/SolverInfo.h>
#include <jive/solver/SolverParams.h>
#include <jive/solver/SolverFactory.h>
#include <jive/solver/SolverException.h>
#include <jive/solver/StdConstrainer.h>
#include <jive/solver/DummyConstrainer.h>

#include "AMGXSolver.h"


JEM_DEFINE_CLASS( jive::solver::AMGXSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::System;
using jem::newInstance;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class AMGXSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   AMGXSolver::TYPE_NAME          = "AmgX";

const int     AMGXSolver::NEW_VALUES_        = 1 << 0;
const int     AMGXSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


AMGXSolver::AMGXSolver

  ( const String&        name,
    Ref<AbstractMatrix>  matrix,
    Ref<Constraints>     cons ) :

    Super ( name )

{
  using jem::util::connect;
  using jive::util::joinNames;

  JEM_PRECHECK ( matrix );

  String  conName = joinNames ( myName_, "constrainer" );


  if ( ! matrix->hasExtension<SparseMatrixExt>() )
  {
    throw jem::IllegalArgumentException (
      JEM_FUNC,
      matrix->getContext() +
      " does not implement the sparse matrix extension"
    );
  }

  if ( matrix->isDistributed() )
  {
    throw jem::IllegalInputException (
      getContext (),
      getContext () + " does not support distributed matrices"
    );
  }

  if ( ! cons )
  {
    conman_ =

      newInstance<DummyConstrainer> ( conName, matrix );
  }
  else
  {
    conman_ =

      newInstance<StdConstrainer>   ( conName, cons, matrix );
  }

  matrix_    = conman_->getOutputMatrix ();
  mode_      = 0;
  precision_ = 1.e-10; //PRECISION;
  iiter_     = 0;
  error_     = 0.0;
  events_    = ~0x0;
  started_   = 0;

  AMGX_SAFE_CALL(AMGX_initialize());
  AMGX_SAFE_CALL(AMGX_install_signal_handler());

  AMGX_SAFE_CALL(AMGX_config_create_from_file(&AMGXcfg_, "./solver.json"));
  AMGX_SAFE_CALL(AMGX_resources_create_simple(&AMGXrsrc_, AMGXcfg_));

  AMGX_SAFE_CALL(AMGX_matrix_create(&AMGXA_, AMGXrsrc_, AMGX_mode_dDDI));
  AMGX_SAFE_CALL(AMGX_vector_create(&AMGXb_, AMGXrsrc_, AMGX_mode_dDDI));
  AMGX_SAFE_CALL(AMGX_vector_create(&AMGXx_, AMGXrsrc_, AMGX_mode_dDDI));
  
  AMGX_SAFE_CALL(AMGX_solver_create(&AMGXsolver_,    AMGXrsrc_, 
                                     AMGX_mode_dDDI, AMGXcfg_ ));

  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


AMGXSolver::~AMGXSolver ()
{
  // Destroy AMGX resources

  AMGX_SAFE_CALL(AMGX_solver_destroy   ( AMGXsolver_));
  AMGX_SAFE_CALL(AMGX_vector_destroy   ( AMGXx_     ));
  AMGX_SAFE_CALL(AMGX_vector_destroy   ( AMGXb_     ));
  AMGX_SAFE_CALL(AMGX_matrix_destroy   ( AMGXA_     ));
  AMGX_SAFE_CALL(AMGX_resources_destroy( AMGXrsrc_  ));

  AMGX_SAFE_CALL(AMGX_config_destroy( AMGXcfg_ ));
  AMGX_SAFE_CALL(AMGX_finalize());
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void AMGXSolver::start ()
{
  conman_->update ();

  if ( ! started_ )
  {
    matrix_->resetEvents ();
  }

  if ( events_ )
  {
    update_ ();
  }

  started_++;
}


//-----------------------------------------------------------------------
//   finish
//-----------------------------------------------------------------------


void AMGXSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }
}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void AMGXSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void AMGXSolver::improve

  ( const Vector&  lhs,
    const Vector&  rhs )

{
  using jem::dot;
  using jem::max;
  using jem::Float;
  using jem::ArithmeticException;

  JEM_PRECHECK ( lhs.size() == rhs.size() );

  SolverScope  scope    ( *this );

  const idx_t  dofCount = matrix_->size(0);

  Vector       du;
  Vector       r;
  Vector       u;
  Vector       f;

  double       rscale;

  if ( lhs.size() != dofCount )
  {
    util::sizeError ( getContext(),
                      "lhs vector", lhs.size(), dofCount );
  }

  if ( mode_ & PRECON_MODE )
  {
    Matrix  vbuf ( dofCount, 2 );

    r .ref ( vbuf[0] );
    du.ref ( vbuf[1] );
    u .ref ( lhs );
    f .ref ( rhs );
  }
  else
  {
    Matrix  vbuf ( dofCount, 4 );

    r .ref ( vbuf[0] );
    du.ref ( vbuf[1] );
    u .ref ( vbuf[2] );
    f .ref ( vbuf[3] );

    conman_->getRhs  ( f, rhs );
    conman_->initLhs ( u, lhs );
  }

  iiter_ = 0;
  error_ = 0.0;
  rscale = std::sqrt ( dot( f, f ) );

  if ( Float::isNaN( rscale ) )
  {
    throw ArithmeticException (
      getContext (),
      "invalid norm of right-hand vector: NaN"
    );
  }

  if ( Float::isTiny( rscale ) )
  {
    u = 0.0;

    if ( ! (mode_ & PRECON_MODE) )
    {
      conman_->getLhs ( lhs, u );
    }

    return;
  }

  rscale = 1.0 / rscale;

  matrix_->matmul ( r, u );

  r = f - r;

  // Pin memory for the residual and the solution vector

  AMGX_SAFE_CALL(AMGX_pin_memory( r. addr(), sizeof(double) * ( r.size()  )));
  AMGX_SAFE_CALL(AMGX_pin_memory( du.addr(), sizeof(double) * ( du.size() )));

  while ( iiter_ < 1 )
  {

    // Upload residual to GPU device

    AMGX_SAFE_CALL(AMGX_vector_upload  ( AMGXb_, r. size(), 1, r.addr() ));
    AMGX_SAFE_CALL(AMGX_vector_set_zero( AMGXx_, du.size(), 1 ));

    // Solve 

    AMGX_SAFE_CALL(AMGX_solver_solve_with_0_initial_guess(AMGXsolver_, AMGXb_, AMGXx_));

    // Download the solution vector

    AMGX_SAFE_CALL(AMGX_vector_download(AMGXx_, du.addr() ));

    u  += du;

    matrix_->matmul ( r, u );

    r = f - r;

    iiter_++;
    error_ = rscale * std::sqrt ( dot( r, r ) );

    solveEvent.emit ( error_, *this );

    if ( Float::isNaN( error_ ) )
    {
      throw ArithmeticException (
        getContext (),
        "invalid norm of residual vector: NaN"
      );
    }

    if ( error_ <= precision_ || error_ > 1.0e5 )
    {
      break;
    }
  }

  if ( (error_ > max( 1.0, precision_ )) ||
       (error_ > precision_ && ! (mode_ & LENIENT_MODE)) )
  {
    throw SolverException (
      getContext     (),
      String::format ( "residual norm too large: %e", error_ )
    );
  }

  if ( ! (mode_ & PRECON_MODE) )
  {
    conman_->getLhs ( lhs, u );
  }

  AMGX_SAFE_CALL(AMGX_unpin_memory( offsets_.addr() ));
  AMGX_SAFE_CALL(AMGX_unpin_memory( indices_.addr() ));
  AMGX_SAFE_CALL(AMGX_unpin_memory( values_ .addr() ));
  
  AMGX_SAFE_CALL(AMGX_unpin_memory( r. addr() ));
  AMGX_SAFE_CALL(AMGX_unpin_memory( du.addr() ));
}


//-----------------------------------------------------------------------
//   getInfo
//-----------------------------------------------------------------------


void AMGXSolver::getInfo ( const Properties& info ) const
{
  double  memUsage = 0.0;
  idx_t   dofCount = matrix_->size (0);

  /*if ( data_ )
  {
    memUsage = data_->getMemUsage ();
  }*/

  Super::getInfo ( info );

  info.set ( SolverInfo::TYPE_NAME,  TYPE_NAME );
  info.set ( SolverInfo::MEM_USAGE,  memUsage  );
  info.set ( SolverInfo::RESIDUAL,   error_    );
  info.set ( SolverInfo::ITER_COUNT, iiter_    );
  info.set ( SolverInfo::DOF_COUNT,  dofCount  );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void AMGXSolver::configure ( const Properties& props )
{}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void AMGXSolver::getConfig ( const Properties& conf ) const
{}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void AMGXSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int AMGXSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void AMGXSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double AMGXSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* AMGXSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* AMGXSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void AMGXSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void AMGXSolver::setZeroThreshold ( double eps )
{
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double AMGXSolver::getZeroThreshold () const
{
  return NAN;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void AMGXSolver::setMaxZeroPivots ( idx_t n )
{
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t AMGXSolver::getMaxZeroPivots () const
{
  return -1;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> AMGXSolver::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  params,
    const Properties&  globdat )

{
  Ref<AbstractMatrix>  mat;
  Ref<Constraints>     cons;


  params.find ( mat,  SolverParams::MATRIX      );
  params.find ( cons, SolverParams::CONSTRAINTS );

  if ( mat && mat->hasExtension<SparseMatrixExt>() )
  {
    return newInstance<Self> ( name, mat, cons );
  }
  else
  {
    return nullptr;
  }
}


//-----------------------------------------------------------------------
//   declare
//-----------------------------------------------------------------------


void AMGXSolver::declare ()
{
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void AMGXSolver::update_ ()
{
  using jem::maxOf;
  using jem::castTo;
  using jem::iarray;
  using jem::System;
  using jem::ArithmeticException;
  using jive::solver::SolverException;

  const String  context    = getContext ();

  SparseMatrixExt*  sx     =

    matrix_->getExtension<SparseMatrixExt> ();

  SparseMatrix      sm     = sx->toSparseMatrix ();

  const idx_t       msize  = sm.size (0);
  const idx_t       nnz    = sm.nonZeroCount ();

  debug_ = & System::debug ( myName_ );

  if ( sm.size(0) != sm.size(1) )
  {
    util::nonSquareMatrixError ( getContext(), sm.shape() );
  }

  JEM_PRECHECK ( msize <= maxOf<int>() );  

  offsets_.resize ( msize + 1 );
  indices_.resize ( nnz );
  values_ .resize ( nnz );

  offsets_ = castTo<int> ( sm.getRowOffsets() );
  indices_ = castTo<int> ( sm.getColumnIndices() );
  values_  =               sm.getValues ();

  // Pin memory

  AMGX_SAFE_CALL(AMGX_pin_memory(offsets_.addr() , sizeof(int)    * (msize + 1) ));
  AMGX_SAFE_CALL(AMGX_pin_memory(indices_.addr() , sizeof(int)    * (nnz) ));
  AMGX_SAFE_CALL(AMGX_pin_memory(values_ .addr() , sizeof(double) * (nnz) ));
  
  // Upload CSR matrix to GPU device

  AMGX_SAFE_CALL(AMGX_matrix_upload_all( AMGXA_, msize, nnz, 1, 1, 
                                         offsets_.addr(), 
                                         indices_.addr(), 
                                         values_.addr(), 
                                         NULL));

  // Setup AMGX solver

  AMGX_SAFE_CALL(AMGX_solver_setup(AMGXsolver_, AMGXA_));

  matrix_->resetEvents ();

  events_ = 0;

  }


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void AMGXSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void AMGXSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void AMGXSolver::setEvents_ ( int events )
{
  if ( started_ )
  {
    throw jem::IllegalOperationException (
      getContext (),
      "matrix changed while solving a linear system of equations"
    );
  }

  events_ |= events;
}


JIVE_END_PACKAGE( solver )

#endif