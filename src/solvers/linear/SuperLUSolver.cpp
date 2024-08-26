
/** @file SuperLUSolver.cpp
 *  @brief Wrapper to solve a linear system with SuperLU_MT.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 26 August 2024
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2024]
 *       
 */

#if defined(WITH_SUPERLUMT)

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
#include <jive/algebra/DiagMatrixExtension.h>
#include <jive/solver/Names.h>
#include <jive/solver/SolverInfo.h>
#include <jive/solver/SolverParams.h>
#include <jive/solver/SolverFactory.h>
#include <jive/solver/SolverException.h>
#include <jive/solver/StdConstrainer.h>
#include <jive/solver/DummyConstrainer.h>

#include <unordered_set>
#include "SuperLUSolver.h"

JEM_DEFINE_CLASS( jive::solver::SuperLUSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::newInstance;
using jem::util::Flex;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class SuperLUSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   SuperLUSolver::TYPE_NAME          = "SuperLU";
const double  SuperLUSolver::PIVOT_THRESHOLD    = 0.1;
const int     SuperLUSolver::MAX_ITER           = 2;
const char*   SuperLUSolver::NUM_THREADS        = "numThreads";

const char*   SuperLUSolver::REORDER_METHODS[3] =
{
  "None",
  "Matrix",
  "Columns"
};

const int     SuperLUSolver::NEW_VALUES_        = 1 << 0;
const int     SuperLUSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


SuperLUSolver::SuperLUSolver

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

  // Set default SuperLU options
  nThreads_  = 1;
  cPermSpec_ = 3;
  panelSize_ = sp_ienv(1);
  relax_     = sp_ienv(2);
  refact_    = NO;
  usePr_     = NO;
  fact_      = EQUILIBRATE;
  trans_     = TRANS;
  diagPivThr_= 1.0;
  dropTol_   = 0.0;
  work_      = NULL;
  lwork_     = 0;
  info_      = 0;

  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


SuperLUSolver::~SuperLUSolver ()
{
  pxgstrf_finalize(&sluOpts_, &AC_);
  Destroy_SuperMatrix_Store(&A_);
  StatFree(&Gstat_);
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void SuperLUSolver::start ()
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


void SuperLUSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }
}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void SuperLUSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void SuperLUSolver::improve

  ( const Vector&  lhs,
    const Vector&  rhs )

{
  using jem::dot;
  using jem::max;
  using jem::Float;
  using jem::System;
  using jem::castTo;
  using jem::ArithmeticException;

  JEM_PRECHECK ( lhs.size() == rhs.size() );

  SolverScope  scope    ( *this );

  const int  dofCount = matrix_->size(0);

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
    // ----------------- SuperLU solver starts here! ----------------------

    // Create a right hand side matrix
    dCreate_Dense_Matrix( &B_, dofCount, 1, r.addr(), dofCount,
                          SLU_DN, SLU_D, SLU_GE);

    // Solve
    info_ = 0;
    dgstrs ( trans_, &L_, &U_, rPerm_.data(), 
             cPerm_.data(), &B_, &Gstat_, &info_ );

    if ( info_ == 0 )
    {
      // Extract the solution vector x
      double *sol = (double*)((DNformat*)B_.Store)->nzval;
      for (int i = 0; i < dofCount; ++i) {
          du[i] = sol[i];
      }
    }
    else
    {
      throw ArithmeticException (
        getContext (),
        String::format ("Failed solve, ERROR CODE: %d", info_ )
        );
      exit(info_);
    }

    // ----------------- SuperLU solver ends here! ----------------------

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


void SuperLUSolver::getInfo ( const Properties& info ) const
{
  double  memUsage = 0.0;
  idx_t   dofCount = matrix_->size (0);

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


void SuperLUSolver::configure ( const Properties& props )
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

    myProps.find ( nThreads_,  NUM_THREADS       );
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void SuperLUSolver::getConfig ( const Properties& conf ) const
{
  using jem::util::setBool;

  Properties  myConf = conf.makeProps ( myName_ );

  Super::getConfig ( conf );

  myConf.set ( PropNames::PIVOT_THRESHOLD, smallPiv_ );

  myConf.set ( PropNames::REORDER,
               REORDER_METHODS[reorder_] );

  setBool    ( myConf,   PropNames::PRINT_PIVOTS,
               options_, PRINT_PIVOTS );

  myConf.set ( NUM_THREADS,        nThreads_  );
}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void SuperLUSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int SuperLUSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void SuperLUSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double SuperLUSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* SuperLUSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* SuperLUSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void SuperLUSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void SuperLUSolver::setZeroThreshold ( double eps )
{
  small_ = eps;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double SuperLUSolver::getZeroThreshold () const
{
  return small_;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void SuperLUSolver::setMaxZeroPivots ( idx_t n )
{
  maxZeroes_ = n;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t SuperLUSolver::getMaxZeroPivots () const
{
  return maxZeroes_;
}


//-----------------------------------------------------------------------
//   setPivotThreshold
//-----------------------------------------------------------------------


void SuperLUSolver::setPivotThreshold ( double alpha )
{
  smallPiv_ = alpha;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getPivotThreshold
//-----------------------------------------------------------------------


double SuperLUSolver::getPivotThreshold () const
{
  return smallPiv_;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


SuperLUSolver::Options SuperLUSolver::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void SuperLUSolver::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void SuperLUSolver::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getReorderMethod
//-----------------------------------------------------------------------


SuperLUSolver::ReorderMethod SuperLUSolver::getReorderMethod () const
{
  return reorder_;
}


//-----------------------------------------------------------------------
//   setReorderMethod
//-----------------------------------------------------------------------


void SuperLUSolver::setReorderMethod ( ReorderMethod method )
{
  reorder_ = method;

  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> SuperLUSolver::makeNew

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


void SuperLUSolver::declare ()
{
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void SuperLUSolver::update_ ()
{
  using jem::isTiny;
  using jem::maxOf;
  using jem::castTo;
  using jem::iarray;
  using jem::System;
  using jem::ArithmeticException;
  using jive::algebra::DiagMatrixExt;
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

  offsets_.resize  ( msize + 1 ); // row offsets
  jindices_.resize ( nnz );       // column indices
  values_ .resize  ( nnz );       // values

  offsets_  = castTo<int> ( sm.getRowOffsets()    );
  jindices_ = castTo<int> ( sm.getColumnIndices() );
  values_   =               sm.getValues ();

  // Set up SuperLU stats
  StatAlloc(msize, nThreads_, panelSize_, relax_, &Gstat_);
  StatInit (msize, nThreads_, &Gstat_);

  // Create a SuperLU matrix of type CSC (NC), although
  // Jive uses a CSR (NR) format.

  dCreate_CompCol_Matrix(&A_, msize, msize, nnz, 
                         values_.addr(),
                         jindices_.addr(), 
                         offsets_.addr(), 
                         SLU_NC, SLU_D, SLU_GE);

  // Get the permutations

  cPerm_.resize(msize);
  rPerm_.resize(msize);

  get_perm_c( cPermSpec_, &A_, cPerm_.data() );

  // Apply permutation to matrix A and store output
  // in matrix AC. Also, initialize the sluOpts.

  pdgstrf_init( nThreads_, fact_, trans_, refact_, panelSize_,
                relax_, diagPivThr_, usePr_, dropTol_, cPerm_.data(),
                rPerm_.data(), work_, lwork_, &A_, &AC_, &sluOpts_, 
                &Gstat_ );

  // Compute the LU factorization of matrix AC
  info_ = 0;
  pdgstrf( &sluOpts_, &AC_, rPerm_.data(), 
           &L_, &U_, 
           &Gstat_, &info_ );

  if ( info_ != 0 )
  {
    throw ArithmeticException (
      getContext (),
      String::format ("Failed factorization, ERROR CODE: %d", info_ )
      );
    exit(info_);
  }

  updated_ = true;  

  matrix_->resetEvents ();

  events_ = 0;

  }


//-----------------------------------------------------------------------
//   connectToSolver_
//-----------------------------------------------------------------------


void SuperLUSolver::connectToSolver_ ()
{
}


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void SuperLUSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void SuperLUSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void SuperLUSolver::setEvents_ ( int events )
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


void SuperLUSolver::progressHandler_ ( idx_t jcol )
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


void SuperLUSolver::pivotHandler_ ( idx_t irow, double pivot )
{
}


//-----------------------------------------------------------------------
//   zeroPivotHandler_
//-----------------------------------------------------------------------


void SuperLUSolver::zeroPivotHandler_ ( idx_t irow, double pivot )
{
}


JIVE_END_PACKAGE( solver )

#endif