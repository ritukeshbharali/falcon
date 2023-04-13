
/** @file PardisoSolver.cpp
 *  @brief Wrapper to solve a linear system with Intel MKL Pardiso.
 *  
 *  This class implements a wrapper to solve the linear
 *  system of equations using the Intel Pardiso solver.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 03 May 2022
 * 
 *  NOTE: Intel oneapi must be installed in the system!
 *        The compiler would give a warning about a
 *        conflict in C language declaration. This is
 *        because void, int are used instead of the mkl
 *        type _MKL_DSS_HANDLE_t, MKL_INT. This should
 *        not be a problem during computation.
 * 
 *        Locate and source "setvars.sh" (typically in
 *        folder /opt/intel/oneapi). This sets up the
 *        necessary environment MKLROOT, prior to 
 *        compiling the executable and running the 
 *        program.
 *
 *        * IMPORTANT * Set sortColumns = 1 for higher
 *        order elements.
 * 
 *
 *  Updates (when, what and who)
 *     - [04 May 2022] Throw error for wrong matrix type,
 *       option to print Pardiso info to terminal with
 *       msglvl set to 0 or 1. (RB)
 * 
 *     - [15 May 2022] Columns must be sorted for each
 *       row in increasing order for Pardiso. For this
 *       code, sorting is required only for higher order
 *       elements. Implemented sortColumns (0 or 1) as
 *       off/on. (RB)
 * 
 *     - [09 Feb 2023] Matrix checker was on by default.
 *       It is now a user-defined option, to be used for
 *       debugging. Changed iparm[1] from Metis Nested
 *       Dissection (2) to Parallel Nested Dissection (3)
 *       It decreases the computational time. (RB)
 *       
 */


/* Include jem and jive headers */

#include <cmath>
#include <vector>
#include <algorithm>
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
#include "PardisoSolver.h"


JEM_DEFINE_CLASS( jive::solver::PardisoSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::newInstance;
using jem::util::Flex;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class PardisoSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   PardisoSolver::TYPE_NAME          = "Pardiso";
const double  PardisoSolver::PIVOT_THRESHOLD    = 0.1;
const int     PardisoSolver::MAX_ITER           = 10;

const char*   PardisoSolver::MATRIX_TYPE        = "mtype";
const char*   PardisoSolver::PRINT_INFO         = "msglvl";
const char*   PardisoSolver::NUM_THREADS        = "numThreads";
const char*   PardisoSolver::SORT_COLUMNS       = "sortColumns";
const char*   PardisoSolver::MATRIX_CHECKER     = "matrixChecker";

const char*   PardisoSolver::REORDER_METHODS[3] =
{
  "None",
  "Matrix",
  "Columns"
};

const int     PardisoSolver::NEW_VALUES_        = 1 << 0;
const int     PardisoSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


PardisoSolver::PardisoSolver

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


  /* Initialize the Pardiso solver  */

  for ( int i = 0; i < 64; i++ )
  {
    iparm_[i] = 0;
  }

  /**
   * List of key input output parameters (iparm_[id])
   * 
   * ------------------------------------------------------------------------------------------------------------------ 
   * id       TYPE       FEATURE                         OPTIONS
   * ------------------------------------------------------------------------------------------------------------------ 
   * 
   * 0        IN         Configuration                   0: default, 1: user
   * 1        IN         Reordering matrix               0:amd, 1:METIS nd, 2: parallel nd
   * 3        IN         Preconditioned CGS/CG           31,61: nonsym, 62: sym
   * 4        IN         Matrix permutation              0:default, 1: user, 2: perm phase 1
   * 5        IN         Solution storage                0: du, 1: rhs replaced by solution
   * 6        OUT        Iterative refinement steps
   * 7        IN         Iterative refinement            0: 2 steps, >0: as specified
   * 8        IN         Relative res. tol.           
   * 9        IN         Pivoting perturbation           (default) 13: unsym, 8: sym
   * 10       IN         Scaling vectors                 0: off (sym indef), 1: on
   * 11       IN         System format Ax=b              0: same, 1: conjugate, 2: transpose
   * 12       IN         Weighted mateching              0: off (sym indef), 1: on
   * 13       OUT        Number of perturbed pivots
   * 14       OUT        Peak memory on symbolic 
   *                     factorization
   * 15       OUT        Permanent memory on symbolic 
   *                     factorization
   * 16       OUT        Size of factors/Peak memory 
   *                     on numerical factorization 
   *                     and solution. 
   * 17       IN/OUT     Report nnz elements             <0: on, >=0 off
   * 18       IN/OUT     Report fp operations in 1e6     <0: on, >=0 off
   * 19       IN/OUT     Report CG diagnostics
   * 20       IN/OUT     Pivoting sym indefinite matrix  0: 1x1 diag, 1,2,3: Bunch-Kaufman
   * 21       OUT        Inertia: number of positive 
   *                     eigenvalues
   * 22       OUT        Inertia: number of negative 
   *                     eigenvalues
   * 23       IN         Parallel factorization control  0: classic, 1: two-level (higher threads > 8), 10: unsym
   * 24       IN         Parallel forward/backward solve 0: matrix or num of rhs partition, 1: sequential, 2: only matrix partition
   * 26       IN         Matrix checker                  0: off, 1: on
   * 27       IN         Precision                       0: double, 1: single
   * 29       OUT        Number of zero or negative 
   *                     pivots
   * 30       IN         Partial solve                   0: off, 1,2,3: on
   * 33       IN         Optimal number of OpenMP threads 
   *                     for conditional numerical 
   *                     reproducibility (CNR) mode
   * 34       IN         Indexing                        0: Fortran indexing, 1: C indexing
   * 35       IN         Schur                           0: off, 1: on
   * 36       IN         Matrix type                     0: CSR matrix, 1: BSR matrix
   * 38       IN         Low rank update                 0: off, 1: on
   * 42       IN         Diagonal of inverse matrix      0: off, 1: on
   * 55       IN         Diagonal and pivoting           0: off, 1: on
   * 59       IN         Core computation                0: in-core, 1: switch, 2: out-of-core
   * 62       OUT        Size of out-of-core memory
   * ------------------------------------------------------------------------------------------------------------------ 
   */

  iparm_[0]   = 1;            // User-defined 
  iparm_[1]   = 3;            // OpenMP parallel nested dissection  
  iparm_[9]   = 13;           // Default, 13 for unsymm, change to 8 for sym indef
  iparm_[10]  = 1;            // Default, 1 for unsymm, change to 0 for sym indef  
  iparm_[11]  = 0;             // Ax = b, no conjugate/transpose
  iparm_[12]  = 1;            // Enable weighted matching
  iparm_[23]  = 10;           // Parallel factorization for unsym
  iparm_[34]  = 1;            // Zero based indexing

  maxfct_     = 1;            // Maximal number of factors in memory
  mnum_       = 1;            // The number of matrix (from 1 to maxfct) to solve 
  mtype_      = 11;           // Real and non-symmetric
  nrhs_       = 1;            // Number of Right-hand side vectors

  msglvl_     = 0;            // Do not print statistical information
  sortCols_   = 0;            // Switch off sorting columns for each row
  matrixChk_  = 0;            // Switch off matrix checker 
  errorFlag_  = 0;            // Initialize error flag to zero
  nThreads_   = 1;            // Default, use one thread


  /* -------------------------------------------------------------------- */
  /* Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */

  for ( int i = 0; i < 64; i++ )
  {
    pt_[i] = 0;
  }

  /* -------------------------------------------------------------------- */
  /* Delete the work array if the matrix is modified                      */
  /* -------------------------------------------------------------------- */

  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );

}


PardisoSolver::~PardisoSolver ()
{
  freePardisoMemory_ ();
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void PardisoSolver::start ()
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


void PardisoSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }

}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void PardisoSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void PardisoSolver::improve

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

  iparm_[7]  = MAX_ITER;

  if ( nThreads_ > 8 )
  {
    iparm_[23] = 1;  // improves scalability
  }

  iparm_[26]  = matrixChk_;   // Matrix checker

  while ( iiter_ < MAX_ITER )
  {
    
    // Pardiso Solver: Solve, iterative refinement

    phase_ = 33;

    pardiso (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
             &n_, values_.addr(), offsets_.addr(), indices_.addr(), 
             &idum_, &nrhs_, iparm_, &msglvl_, r.addr(), du.addr(), &errorFlag_);

    if ( errorFlag_ != 0 )
    {
    switch(errorFlag_)
    {
      case -1 :
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " input inconsistent!!!"
        );
        break;
      case -2 :
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " not enough memory!!!"
        );
        break;
      case -3 : 
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " reordering problem!!!"
        );
      case -4 : 
        throw jem::ArithmeticException (
          getContext (),
          getContext () + " zero pivot!!!"
        ); 
      case -9 : 
        throw jem::ArithmeticException (
          getContext (),
          getContext () + " insufficient memory for out-of-core!!!"
        );    
        break;
      default :
        throw ArithmeticException (
        getContext (),
        String::format ("Something went wrong, ERROR CODE: %d", errorFlag_ )
      );
        break;
    }
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


void PardisoSolver::getInfo ( const Properties& info ) const
{

  double  memUsage = 0.0;
  // idx_t   dofCount = matrix_->size (0);

  Super::getInfo ( info );

  info.set ( SolverInfo::TYPE_NAME,  TYPE_NAME );
  info.set ( SolverInfo::MEM_USAGE,  memUsage  );
  info.set ( SolverInfo::RESIDUAL,   error_    );
  info.set ( SolverInfo::ITER_COUNT, iiter_    );
  info.set ( SolverInfo::DOF_COUNT,  n_        );

}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void PardisoSolver::configure ( const Properties& props )
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

    myProps.find ( mtype_,     MATRIX_TYPE     );
    myProps.find ( nThreads_,  NUM_THREADS     );
    myProps.find ( msglvl_,    PRINT_INFO      );
    myProps.find ( sortCols_,  SORT_COLUMNS    );
    myProps.find ( matrixChk_, MATRIX_CHECKER  );

    using std::vector;

    vector<int> mtypes = {1,2,-2,3,4,-4,6,11,13};
    vector<int>::iterator it;
    it = std::find (mtypes.begin(), mtypes.end(), mtype_);

    if (it == mtypes.end())
    {
      throw jem::IllegalInputException (
      getContext (),
      getContext () + " Wrong matrix type: Choose 1,2,-2,3,4,-4,6,11,13!"
    );
        
    }
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void PardisoSolver::getConfig ( const Properties& conf ) const
{

  using jem::util::setBool;

  Properties  myConf = conf.makeProps ( myName_ );

  Super::getConfig ( conf );

  myConf.set ( PropNames::PIVOT_THRESHOLD, smallPiv_ );

  myConf.set ( PropNames::REORDER,
               REORDER_METHODS[reorder_] );

  setBool    ( myConf,   PropNames::PRINT_PIVOTS,
               options_, PRINT_PIVOTS );

  myConf.set ( MATRIX_TYPE,    mtype_     );
  myConf.set ( NUM_THREADS,    nThreads_  );
  myConf.set ( PRINT_INFO,     msglvl_    );
  myConf.set ( SORT_COLUMNS,   sortCols_  );
  myConf.set ( MATRIX_CHECKER, matrixChk_ );

  // Switch off Intel dynamic parallelism and set threads

  mkl_set_dynamic(0);
  mkl_set_num_threads(nThreads_);

}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void PardisoSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int PardisoSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void PardisoSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double PardisoSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* PardisoSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* PardisoSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void PardisoSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void PardisoSolver::setZeroThreshold ( double eps )
{
  small_ = eps;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double PardisoSolver::getZeroThreshold () const
{
  return small_;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void PardisoSolver::setMaxZeroPivots ( idx_t n )
{
  maxZeroes_ = n;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t PardisoSolver::getMaxZeroPivots () const
{
  return maxZeroes_;
}


//-----------------------------------------------------------------------
//   setPivotThreshold
//-----------------------------------------------------------------------


void PardisoSolver::setPivotThreshold ( double alpha )
{
  smallPiv_ = alpha;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getPivotThreshold
//-----------------------------------------------------------------------


double PardisoSolver::getPivotThreshold () const
{
  return smallPiv_;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


PardisoSolver::Options PardisoSolver::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void PardisoSolver::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void PardisoSolver::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getReorderMethod
//-----------------------------------------------------------------------


PardisoSolver::ReorderMethod PardisoSolver::getReorderMethod () const
{
  return reorder_;
}


//-----------------------------------------------------------------------
//   setReorderMethod
//-----------------------------------------------------------------------


void PardisoSolver::setReorderMethod ( ReorderMethod method )
{
  reorder_ = method;

  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> PardisoSolver::makeNew

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


void PardisoSolver::declare ()
{
  
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void PardisoSolver::update_ ()
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
  n_                       = msize;

  debug_ = & System::debug ( myName_ );

  if ( sm.size(0) != sm.size(1) )
  {
    util::nonSquareMatrixError ( getContext(), sm.shape() );
  }

  JEM_PRECHECK ( msize <= maxOf<int>() );  

  offsets_.resize ( msize + 1 );  // ia
  indices_.resize ( nnz );        // ja
  values_ .resize ( nnz );        // a

  offsets_ = castTo<int> ( sm.getRowOffsets() );
  indices_ = castTo<int> ( sm.getColumnIndices() );
  values_  =               sm.getValues ();

  /* Indices in each row sorted in increasing order */

  if ( sortCols_ == 1 )
  {
    Vector tmp ( n_ );

    for ( idx_t irow = 0; irow < n_; irow++ )
    {
      int  i = offsets_[irow];
      int  k = offsets_[irow + 1];

      for ( int j = i; j < k; j++ )
      {
        tmp[indices_[j]] = values_[j];
      }

      sort ( indices_[slice(i,k)] );

      for ( int j = i; j < k; j++ )
      {
        values_[j] = tmp[indices_[j]];
      }
    }

  }

  /* Symbolic factorization */ 

  phase_ = 11;

  pardiso (  pt_, &maxfct_, &mnum_, &mtype_, &phase_,
             &n_, values_.addr(), offsets_.addr(), indices_.addr(), 
             &idum_, &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &errorFlag_ );

  if ( errorFlag_ != 0 )
  {
    throw ArithmeticException (
      getContext (),
      String::format ("Symbolic factorization, ERROR CODE: %d", errorFlag_ )
      );

  }

  /* Numeric factorization */ 

  phase_ = 22;

  pardiso (  pt_, &maxfct_, &mnum_, &mtype_, &phase_,
             &n_, values_.addr(), offsets_.addr(), indices_.addr(), 
             &idum_, &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &errorFlag_ );

  if ( errorFlag_ != 0 )
  {
    throw ArithmeticException (
      getContext (),
      String::format ("Numeric factorization, ERROR CODE: %d", errorFlag_ )
      );

  }

  matrix_->resetEvents ();

  events_ = 0;

}


//-----------------------------------------------------------------------
//   connectToSolver_
//-----------------------------------------------------------------------


void PardisoSolver::connectToSolver_ ()
{
}


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void PardisoSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void PardisoSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void PardisoSolver::setEvents_ ( int events )
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
//   freePardisoMemory_
//-----------------------------------------------------------------------


void PardisoSolver::freePardisoMemory_ ()
{

  phase_ = -1;

  pardiso (  pt_, &maxfct_, &mnum_, &mtype_, &phase_,
             &n_, values_.addr(), offsets_.addr(), indices_.addr(), 
             &idum_, &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &errorFlag_ );



}


//-----------------------------------------------------------------------
//   progressHandler_
//-----------------------------------------------------------------------


void PardisoSolver::progressHandler_ ( idx_t jcol )
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


void PardisoSolver::pivotHandler_ ( idx_t irow, double pivot )
{
}


//-----------------------------------------------------------------------
//   zeroPivotHandler_
//-----------------------------------------------------------------------


void PardisoSolver::zeroPivotHandler_ ( idx_t irow, double pivot )
{
}


JIVE_END_PACKAGE( solver )
