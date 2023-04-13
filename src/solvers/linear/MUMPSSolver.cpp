
/** @file MUMPSSolver.cpp
 *  @brief Wrapper to solve a linear system with MUMPS (sequential only).
 *  
 *  This class implements a wrapper to solve the
 *  linear system of equations using the sequential 
 *  multi-threaded version of MUMPS.
 *  (http://mumps.enseeiht.fr/).
 *  
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 13 May 2022
 * 
 *  NOTE: Requires a sequential MUMPS version. Works
 *        only for the sequential version of jiveFEA.
 * 
 *  TO-DO: Add more user-defined controls. 
 *         Extend to MPI+Threading.
 *
 *  Updates (when, what and who)
 *     - [16 May 2022] Combine Factorization, Solve,
 *       and Back substitution into a single call to
 *       MUMPS, instead of three separate calls. (RB)
 *       
 */


/* Include jem and jive headers */

#include <cmath>
#include <jem/base/assert.h>
#include <jem/base/limits.h>
#include <jem/base/Array.h>
#include <jem/base/Float.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/ArithmeticException.h>
#include <jem/base/IllegalInputException.h>
#include <jem/base/IllegalArgumentException.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/base/array/operators.h>
#include <jem/base/array/utilities.h>
#include <jem/util/Flex.h>
#include <jem/numeric/sparse/Reorder.h>
#include <jem/numeric/sparse/SparseLU.h>
#include <jem/util/Event.h>
#include <jive/util/error.h>
#include <jive/util/utilities.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/SparseMatrixObject.h>
#include <jive/solver/Names.h>
#include <jive/solver/SolverInfo.h>
#include <jive/solver/SolverParams.h>
#include <jive/solver/SolverFactory.h>
#include <jive/solver/SolverException.h>
#include <jive/solver/StdConstrainer.h>
#include <jive/solver/DummyConstrainer.h>

/* MUMPS */

#include "MUMPSSolver.h"

#define JOB_INIT   -1
#define JOB_END    -2
#define JOB_SOLVE   3
#define JOB_FACTOR  4
#define JOB_ALL     6
#define USE_COMM_WORLD -987654

#define ICNTL(I) icntl[(I)-1]   /* macro s.t. indices match documentation */


JEM_DEFINE_CLASS( jive::solver::MUMPSSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::newInstance;
using jem::util::Flex;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class MUMPSSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   MUMPSSolver::TYPE_NAME          = "MUMPS";
const double  MUMPSSolver::PIVOT_THRESHOLD    = 0.1;
const int     MUMPSSolver::MAX_ITER           = 1;

const char*   MUMPSSolver::REORDER_METHODS[3] =
{
  "None",
  "Matrix",
  "Columns"
};

const int     MUMPSSolver::NEW_VALUES_        = 1 << 0;
const int     MUMPSSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


MUMPSSolver::MUMPSSolver

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
  small_     = ZERO_THRESHOLD;
  smallPiv_  = PIVOT_THRESHOLD;
  precision_ = PRECISION;
  reorder_   = REORDER_MATRIX;
  options_   = 0;
  maxZeroes_ = 0;
  iiter_     = 0;
  error_     = 0.0;
  events_    = ~0x0;
  started_   = 0;

  /* Initialize dummy MPI from MUMPS (libmpiseq) */

  int argc = 1;
  char ** argv;

  mpierr_ = MPI_Init( &argc, &argv );
  mpierr_ = MPI_Comm_rank( MPI_COMM_WORLD, &mpirank_ );

  /* Setup MUMPS parameters  */

  id_.par = 1;                        // Host performs factorization, solve
  id_.sym = 0;                        // Unsymmetric matrix
  id_.comm_fortran = USE_COMM_WORLD;  // Default (as per manual)

  /* Initialize a MUMPS Instance */

  id_.job = JOB_INIT;
  dmumps_c ( &id_ );

  /* Check for errors */

  if ( id_.infog[0] < 0 )
  {
    throw SolverException (
      getContext     (),
      String::format ( "MUMPS Error: INFOG(1) : %d, INFOG(2): %d", id_.infog[0], id_.infog[1] )
      );
  }
  
  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


MUMPSSolver::~MUMPSSolver ()
{
  id_.job = JOB_END;
  mpierr_ = MPI_Finalize();
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void MUMPSSolver::start ()
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


void MUMPSSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }
}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void MUMPSSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void MUMPSSolver::improve

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

  while ( iiter_ < MAX_ITER )
  {
    
    // MUMPS replaces the residual with solution
    // So, residual is copied to solution vector
    // to keep the jive structure for the solver.

    du = r;

    // Solve with MUMPS

    id_.job  = JOB_ALL;
    id_.rhs  = du.addr();
    id_.nrhs = 1;
    id_.lrhs = id_.n;

    dmumps_c( &id_ );

    /* Check for errors */

    if ( id_.infog[0] < 0 )
    {
      throw SolverException (
        getContext     (),
        String::format ( "MUMPS Error: INFOG(1) : %d, INFOG(2): %d", id_.infog[0], id_.infog[1] )
        );
    }

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
}


//-----------------------------------------------------------------------
//   getInfo
//-----------------------------------------------------------------------


void MUMPSSolver::getInfo ( const Properties& info ) const
{
  double  memUsage = 0.0;
  idx_t   dofCount = matrix_->size (0);

/*  if ( data_ )
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


void MUMPSSolver::configure ( const Properties& props )
{
  using jem::util::findBool;

  Super::configure ( props );

  if ( props.contains( myName_ ) )
  {
    Properties  myProps = props.findProps ( myName_ );

    String      reord;

    myProps.find ( smallPiv_, PropNames::PIVOT_THRESHOLD,
                   0.0, 1.0 );

    if ( myProps.find( reord, PropNames::REORDER ) )
    {
      int  i;

      for ( i = 0; i < 3; i++ )
      {
        if ( reord.equalsIgnoreCase( REORDER_METHODS[i] ) )
        {
          break;
        }
      }

      if ( i >= 3 )
      {
        myProps.propertyError (
          PropNames::REORDER,
          String::format (
            "invalid re-ordering method; "
            "should be `%s\', `%s\' or `%s\'",
            REORDER_METHODS[0],
            REORDER_METHODS[1],
            REORDER_METHODS[2]
          )
        );
      }

      reorder_ = (ReorderMethod) i;
    }

    findBool ( options_, PRINT_PIVOTS,
               myProps,  PropNames::PRINT_PIVOTS );

  }

  id_.ICNTL(1)  = -1;    // suppress error messages
  id_.ICNTL(2)  = -1;    // suppress diagnostic
  id_.ICNTL(3)  = -1;    // suppress global info
  id_.ICNTL(4)  = -1;    // no output messages
  id_.ICNTL(5)  =  0;    // assembled matrix
  id_.ICNTL(6)  =  7;    // automatic zero-free diagonal
  id_.ICNTL(7)  =  7;    // automatic permutation
  id_.ICNTL(8)  =  77;   // automatic scaling
  id_.ICNTL(9)  =  1;    // Ax = b
  id_.ICNTL(10) =  2;    // Iterative refinement steps
  id_.ICNTL(11) =  0;    // no error analysis (for performance)
  id_.ICNTL(12) =  0;    // automatic ordering for symmetric matrix
  id_.ICNTL(13) =  0;    // parallel factorization on root node
  id_.ICNTL(16) =  0;    // openMP threads based on application

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MUMPSSolver::getConfig ( const Properties& conf ) const
{
  using jem::util::setBool;

  Properties  myConf = conf.makeProps ( myName_ );

  Super::getConfig ( conf );

  myConf.set ( PropNames::PIVOT_THRESHOLD, smallPiv_ );

  myConf.set ( PropNames::REORDER,
               REORDER_METHODS[reorder_] );

  setBool    ( myConf,   PropNames::PRINT_PIVOTS,
               options_, PRINT_PIVOTS );
}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void MUMPSSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int MUMPSSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void MUMPSSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double MUMPSSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* MUMPSSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* MUMPSSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void MUMPSSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void MUMPSSolver::setZeroThreshold ( double eps )
{
  small_ = eps;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double MUMPSSolver::getZeroThreshold () const
{
  return small_;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void MUMPSSolver::setMaxZeroPivots ( idx_t n )
{
  maxZeroes_ = n;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t MUMPSSolver::getMaxZeroPivots () const
{
  return maxZeroes_;
}


//-----------------------------------------------------------------------
//   setPivotThreshold
//-----------------------------------------------------------------------


void MUMPSSolver::setPivotThreshold ( double alpha )
{
  smallPiv_ = alpha;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getPivotThreshold
//-----------------------------------------------------------------------


double MUMPSSolver::getPivotThreshold () const
{
  return smallPiv_;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


MUMPSSolver::Options MUMPSSolver::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void MUMPSSolver::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void MUMPSSolver::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getReorderMethod
//-----------------------------------------------------------------------


MUMPSSolver::ReorderMethod MUMPSSolver::getReorderMethod () const
{
  return reorder_;
}


//-----------------------------------------------------------------------
//   setReorderMethod
//-----------------------------------------------------------------------


void MUMPSSolver::setReorderMethod ( ReorderMethod method )
{
  reorder_ = method;

  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> MUMPSSolver::makeNew

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


void MUMPSSolver::declare ()
{
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void MUMPSSolver::update_ ()
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

  id_.n   = msize;
  id_.nnz = nnz; 

  debug_ = & System::debug ( myName_ );

  if ( sm.size(0) != sm.size(1) )
  {
    util::nonSquareMatrixError ( getContext(), sm.shape() );
  }

  JEM_PRECHECK ( msize <= maxOf<int>() );  

  offsets_.resize  ( msize + 1 ); // row offsets
  jindices_.resize ( nnz );       // column indices
  values_ .resize  ( nnz );       // values


  offsets_  = castTo<int> ( sm.getRowOffsets()    );
  jindices_ = castTo<int> ( sm.getColumnIndices() );
  values_   =               sm.getValues ();

  // CSR to COO: store row coordinates in offsets_ (Fortran indexing)

  {
    int count = 0;

    jem::Array<int> tmp ( msize + 1 );
    tmp = offsets_;

    offsets_.resize( nnz );
    offsets_ = 0;

    for ( idx_t irow = 0; irow < msize; irow++ )
    {
      for ( idx_t k = 0; k < (tmp[irow+1] - tmp[irow]); k++ )
      {
        offsets_[count] = irow + 1;
        count++;
      }
    }

  }

  // Update column indices to Fortran indexing

  for ( idx_t i = 0; i < nnz; i++ ) jindices_[i]++;

  // Setup COO for MUMPS  

  id_.irn = offsets_. addr();
  id_.jcn = jindices_.addr();
  id_.a   = values_  .addr();

  matrix_->resetEvents ();

  events_ = 0;

  }


//-----------------------------------------------------------------------
//   connectToSolver_
//-----------------------------------------------------------------------


void MUMPSSolver::connectToSolver_ ()
{
}


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void MUMPSSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void MUMPSSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void MUMPSSolver::setEvents_ ( int events )
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


//-----------------------------------------------------------------------
//   progressHandler_
//-----------------------------------------------------------------------


void MUMPSSolver::progressHandler_ ( idx_t jcol )
{
  const idx_t  msize = matrix_->size(0);

  if ( msize <= 1 )
  {
    factorEvent.emit ( 1, * this );
  }
  else
  {
    idx_t  done = (idx_t) ((100.0 * (double) jcol) /
                                    (double) (msize - 1));

    factorEvent.emit ( done, *this );
  }
}


//-----------------------------------------------------------------------
//   pivotHandler_
//-----------------------------------------------------------------------


void MUMPSSolver::pivotHandler_ ( idx_t irow, double pivot )
{
}


//-----------------------------------------------------------------------
//   zeroPivotHandler_
//-----------------------------------------------------------------------


void MUMPSSolver::zeroPivotHandler_ ( idx_t irow, double pivot )
{
}


JIVE_END_PACKAGE( solver )
