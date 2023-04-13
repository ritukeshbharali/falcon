
/** @file UmfpackSolver.cpp
 *  @brief Wrapper to solve a linear system with Umfpack.
 *  
 *  This class implements a wrapper to solve the linear
 *  system of equations using the Umfpack solver. (doi:
 *  10.1145/992200.992206) It is inspired from mtxuss :
 *  the MTX unsymmetric sparse matrix program, written
 *  by Frank Everdij, TU Delft, the Netherlands. 
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 07 April 2022
 * 
 *  NOTE: Umfpack must be installed in the system!
 *  
 *
 *  Updates (when, what and who)
 *     - [XX YYYYY 2022]
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
#include "UmfpackSolver.h"


JEM_DEFINE_CLASS( jive::solver::UmfpackSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::newInstance;
using jem::util::Flex;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class UmfpackSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   UmfpackSolver::TYPE_NAME          = "Umfpack";
const double  UmfpackSolver::PIVOT_THRESHOLD    = 0.1;
const int     UmfpackSolver::MAX_ITER           = 10;

const char*   UmfpackSolver::REORDER_METHODS[3] =
{
  "None",
  "Matrix",
  "Columns"
};

const int     UmfpackSolver::NEW_VALUES_        = 1 << 0;
const int     UmfpackSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


UmfpackSolver::UmfpackSolver

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

  numeric_   = NULL;

  /* Set default Umfpack parameters  */

  umfpack_dl_defaults( control_ );

  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


UmfpackSolver::~UmfpackSolver ()
{
  freeNumeric_ ();
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void UmfpackSolver::start ()
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


void UmfpackSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }
}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void UmfpackSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void UmfpackSolver::improve

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
  int          umferr;

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
    
    // Solve using Umfpack

    umferr = umfpack_dl_solve ( UMFPACK_At,
                                offsets_.addr (),
                                indices_.addr (),
                                values_ .addr (),
                                du      .addr (),
                                r       .addr (),
                                numeric_,
                                control_,
                                info_ );
    if ( umferr < 0 )
    {
      umfpack_dl_report_info(control_, info_);
      umfpack_dl_report_status(control_, umferr);

      throw ArithmeticException (
        getContext (),
        String::format ("error in umfpack_di_solve: %d", umferr )
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


void UmfpackSolver::getInfo ( const Properties& info ) const
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


void UmfpackSolver::configure ( const Properties& props )
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
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void UmfpackSolver::getConfig ( const Properties& conf ) const
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


void UmfpackSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int UmfpackSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void UmfpackSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double UmfpackSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* UmfpackSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* UmfpackSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void UmfpackSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void UmfpackSolver::setZeroThreshold ( double eps )
{
  small_ = eps;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double UmfpackSolver::getZeroThreshold () const
{
  return small_;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void UmfpackSolver::setMaxZeroPivots ( idx_t n )
{
  maxZeroes_ = n;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t UmfpackSolver::getMaxZeroPivots () const
{
  return maxZeroes_;
}


//-----------------------------------------------------------------------
//   setPivotThreshold
//-----------------------------------------------------------------------


void UmfpackSolver::setPivotThreshold ( double alpha )
{
  smallPiv_ = alpha;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getPivotThreshold
//-----------------------------------------------------------------------


double UmfpackSolver::getPivotThreshold () const
{
  return smallPiv_;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


UmfpackSolver::Options UmfpackSolver::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void UmfpackSolver::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void UmfpackSolver::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getReorderMethod
//-----------------------------------------------------------------------


UmfpackSolver::ReorderMethod UmfpackSolver::getReorderMethod () const
{
  return reorder_;
}


//-----------------------------------------------------------------------
//   setReorderMethod
//-----------------------------------------------------------------------


void UmfpackSolver::setReorderMethod ( ReorderMethod method )
{
  reorder_ = method;

  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> UmfpackSolver::makeNew

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


void UmfpackSolver::declare ()
{
  // SolverFactory::declare ( "Umfpack",  & makeNew );
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void UmfpackSolver::update_ ()
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

  int               umferr;
  Vector            accu ( msize );

  debug_ = & System::debug ( myName_ );

  if ( sm.size(0) != sm.size(1) )
  {
    util::nonSquareMatrixError ( getContext(), sm.shape() );
  }

  JEM_PRECHECK ( msize <= maxOf<int>() );  

  offsets_.resize ( msize + 1 );  // Ap
  indices_.resize ( nnz );        // Ai
  values_ .resize ( nnz );        // Ax

  offsets_ = castTo<int> ( sm.getRowOffsets() );
  indices_ = castTo<int> ( sm.getColumnIndices() );
  values_  =               sm.getValues ();

  if ( ! numeric_ )
  {

    /** Jive adopts the Yale sparse format or CSR,
     *  while Umfpack uses a column index format.
     *  So, the matrix entries must be rearranged
     *  according to the column index.
     */ 

    for ( idx_t irow = 0; irow < msize; irow++ )
    {
      int  i = offsets_[irow];
      int  k = offsets_[irow + 1];

      for ( int j = i; j < k; j++ )
      {
        int l   = indices_[j];
        accu[l] = values_[j];
      }

      sort ( indices_[slice(i,k)] );

      for ( int j = i; j < k; j++ )
      {
        int l      = indices_[j];
        values_[j] = accu[l];
      }
    }

    umferr = umfpack_dl_symbolic ( (int) msize, (int) msize,
                                   offsets_.addr (),
                                   indices_.addr (),
                                   values_ .addr (),
                                   &symbolic_,
                                   control_,
                                   info_ );

    if ( umferr < 0 )
    {
      umfpack_dl_report_info(control_, info_);
      umfpack_dl_report_status(control_, umferr);

     throw ArithmeticException (
       getContext (),
       String::format ("error in umfpack_di_symbolic: %d", umferr )
     );
    }

    umferr = umfpack_dl_numeric ( offsets_.addr (),
                                  indices_.addr (),
                                  values_ .addr (),
                                  symbolic_,
                                  &numeric_,
                                  control_,
                                  info_ );

    if ( umferr < 0 )
    {
      umfpack_dl_report_info(control_, info_);
      umfpack_dl_report_status(control_, umferr);

     throw ArithmeticException (
       getContext (),
       String::format ("error in umfpack_di_numeric: %d", umferr )
     );
    }

    freeSymbolic_();

  }

  matrix_->resetEvents ();

  events_ = 0;

  }


//-----------------------------------------------------------------------
//   connectToSolver_
//-----------------------------------------------------------------------


void UmfpackSolver::connectToSolver_ ()
{
}


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void UmfpackSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );

  freeNumeric_();
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void UmfpackSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );

  freeNumeric_();
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void UmfpackSolver::setEvents_ ( int events )
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
//   freeSymbolic_
//-----------------------------------------------------------------------


void UmfpackSolver::freeSymbolic_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack freeSymbolic\n";

  if ( symbolic_ )
  {
    umfpack_dl_free_symbolic ( &symbolic_ );

    symbolic_ = NULL;
  }
}


//-----------------------------------------------------------------------
//   freeNumeric_
//-----------------------------------------------------------------------


void UmfpackSolver::freeNumeric_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack freeNumeric\n";

  if ( numeric_ )
  {
    umfpack_dl_free_numeric ( &numeric_ );

    numeric_ = NULL;
  }
}


//-----------------------------------------------------------------------
//   progressHandler_
//-----------------------------------------------------------------------


void UmfpackSolver::progressHandler_ ( idx_t jcol )
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


void UmfpackSolver::pivotHandler_ ( idx_t irow, double pivot )
{
}


//-----------------------------------------------------------------------
//   zeroPivotHandler_
//-----------------------------------------------------------------------


void UmfpackSolver::zeroPivotHandler_ ( idx_t irow, double pivot )
{
}


JIVE_END_PACKAGE( solver )
