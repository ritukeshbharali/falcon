
/** @file AMGCLSolver.cpp
 *  @brief Wrapper to solve a linear system with AMGCL (sequential only).
 * 
 *  Copyright (C) 2022 Chalmers. All rights reserved.
 *  
 *  This class implements a wrapper to solve the
 *  linear system of equations using the sequential 
 *  multi-threaded version of AMGCL library solvers
 *  and preconditioners.
 *  (https://amgcl.readthedocs.io/en/latest/index.html).
 *  
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 30 June 2023
 * 
 *  @note: Changing the solver and preconditioner based
 *         on user input is not yet implemented. At the
 *         moment, the source file needs to be modified.
 * 
 *  TO-DO: Add more user-defined controls. 
 *         Extend to MPI+Threading.
 *
 *  Updates (when, what and who)
 *       
 */

#if defined(WITH_AMGCL)

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

/* Include AMGCL headers */

#include "AMGCLSolver.h"

#include "../amgcl/amg.hpp"
#include "../amgcl/make_solver.hpp"
#include "../amgcl/adapter/crs_tuple.hpp"
#include "../amgcl/backend/builtin.hpp"
#include "../amgcl/coarsening/aggregation.hpp"
#include "../amgcl/coarsening/plain_aggregates.hpp"
#include "../amgcl/coarsening/pointwise_aggregates.hpp"
#include "../amgcl/coarsening/ruge_stuben.hpp"
#include "../amgcl/coarsening/runtime.hpp"
#include "../amgcl/coarsening/smoothed_aggregation.hpp"
#include "../amgcl/preconditioner/runtime.hpp"
#include "../amgcl/relaxation/as_preconditioner.hpp"
#include "../amgcl/relaxation/damped_jacobi.hpp"
#include "../amgcl/relaxation/gauss_seidel.hpp"
#include "../amgcl/relaxation/ilu0.hpp"
#include "../amgcl/relaxation/iluk.hpp"
#include "../amgcl/relaxation/ilup.hpp"
#include "../amgcl/relaxation/ilut.hpp"
#include "../amgcl/relaxation/spai0.hpp"
#include "../amgcl/relaxation/spai1.hpp"
#include "../amgcl/relaxation/runtime.hpp"
#include "../amgcl/solver/bicgstab.hpp"
#include "../amgcl/solver/cg.hpp"
#include "../amgcl/solver/fgmres.hpp"
#include "../amgcl/solver/gmres.hpp"
#include "../amgcl/solver/idrs.hpp"
#include "../amgcl/solver/lgmres.hpp"
#include "../amgcl/solver/preonly.hpp"
#include "../amgcl/solver/runtime.hpp"

#include "../amgcl/io/mm.hpp"
#include "../amgcl/profiler.hpp"


JEM_DEFINE_CLASS( jive::solver::AMGCLSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::newInstance;
using jem::util::Flex;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class AMGCLSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   AMGCLSolver::TYPE_NAME          = "AMGCL";
const double  AMGCLSolver::PIVOT_THRESHOLD    = 0.1;
const int     AMGCLSolver::MAX_ITER           = 1;
const char*   AMGCLSolver::REORDER_METHODS[3] =
{
  "None",
  "Matrix",
  "Columns"
};

const int     AMGCLSolver::NEW_VALUES_        = 1 << 0;
const int     AMGCLSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


AMGCLSolver::AMGCLSolver

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
  
  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


AMGCLSolver::~AMGCLSolver ()
{
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void AMGCLSolver::start ()
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


void AMGCLSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }
}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void AMGCLSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void AMGCLSolver::improve

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

    // ----------------- AMGCL solver starts here! ----------------------

    // Define AMGCL backend typedefs (required by the library)

    typedef amgcl::backend::builtin<double> SolverBackend;
    typedef amgcl::backend::builtin<double> PreconBackend;

    // Define AMGCL solver and preconditioner

    typedef amgcl::make_solver<
        // Use AMG as preconditioner:
        amgcl::amg<
        PreconBackend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::ilu0
        >,

       // And BiCGStab as iterative solver:
       amgcl::solver::bicgstab<SolverBackend>
       > Solver;

    // Define additional solver parameters
    
    Solver::params prm;
    prm.solver.tol = 1.e-8;   
    prm.precond.coarsening.aggr.eps_strong = 1e-3;

    // Initialize std::vector for solution (x) and 
    // residual (b)

    std::vector<double> x  ( dofCount, 0.0 );
    std::vector<double> b  ( dofCount );

    // Copy the jive Array residual to std::vector b
    #pragma omp parallel for
    for ( int i = 0; i < dofCount; i++ )
    {
      b[i] = r[i];
    }

    // Apply diagonal scaling to the matrix and the residual

    for ( idx_t i = 0; i < dofCount; i++ )
    {
      for ( idx_t j = offsets_[i], e = offsets_[i+1]; j < e; j++ )
      {
        values_[j] *= dscale_[i] * dscale_[jindices_[j]];
      }

      b[i] *= dscale_[i];
    }

    // Tie the CSR Tuple (required by AMGCL)

    auto A = std::tie( dofCount, offsets_, jindices_, values_ );
    
    // Setup the AMGCL solver

    Solver solve( A, prm );

    // Show the mini-report on the constructed solver

    // std::cout << solve << std::endl;

    // Solve the system of equations

    solve(A, b, x);

    // Copy the solution (std::vector) to jive Array
    #pragma omp parallel for
    for ( int i = 0; i < dofCount; i++ )
    {
      du[i] = x[i] * dscale_[i];
    }

    // ----------------- AMGCL solver ends here! ----------------------

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

  /*if ( (error_ > max( 1.0, precision_ )) ||
       (error_ > precision_ && ! (mode_ & LENIENT_MODE)) )
  {
    throw SolverException (
      getContext     (),
      String::format ( "residual norm too large: %e", error_ )
    );
  }*/

  if ( ! (mode_ & PRECON_MODE) )
  {
    conman_->getLhs ( lhs, u );
  }

}


//-----------------------------------------------------------------------
//   getInfo
//-----------------------------------------------------------------------


void AMGCLSolver::getInfo ( const Properties& info ) const
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


void AMGCLSolver::configure ( const Properties& props )
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


void AMGCLSolver::getConfig ( const Properties& conf ) const
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


void AMGCLSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int AMGCLSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void AMGCLSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double AMGCLSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* AMGCLSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* AMGCLSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void AMGCLSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void AMGCLSolver::setZeroThreshold ( double eps )
{
  small_ = eps;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double AMGCLSolver::getZeroThreshold () const
{
  return small_;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void AMGCLSolver::setMaxZeroPivots ( idx_t n )
{
  maxZeroes_ = n;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t AMGCLSolver::getMaxZeroPivots () const
{
  return maxZeroes_;
}


//-----------------------------------------------------------------------
//   setPivotThreshold
//-----------------------------------------------------------------------


void AMGCLSolver::setPivotThreshold ( double alpha )
{
  smallPiv_ = alpha;

  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   getPivotThreshold
//-----------------------------------------------------------------------


double AMGCLSolver::getPivotThreshold () const
{
  return smallPiv_;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


AMGCLSolver::Options AMGCLSolver::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void AMGCLSolver::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void AMGCLSolver::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getReorderMethod
//-----------------------------------------------------------------------


AMGCLSolver::ReorderMethod AMGCLSolver::getReorderMethod () const
{
  return reorder_;
}


//-----------------------------------------------------------------------
//   setReorderMethod
//-----------------------------------------------------------------------


void AMGCLSolver::setReorderMethod ( ReorderMethod method )
{
  reorder_ = method;

  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> AMGCLSolver::makeNew

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


void AMGCLSolver::declare ()
{
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void AMGCLSolver::update_ ()
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

  // Extract diagonal entries of matrix for scaling

  DiagMatrixExt*    dx     = 

    matrix_->getExtension<DiagMatrixExt> ();  

  dscale_.resize   ( msize   );     // diagonals for scaling
  
  dx->getDiagonal  ( dscale_ );

  for ( idx_t i = 0; i < msize; i++ )
  {
    if ( isTiny( dscale_[i] ) )
    {
      dscale_[i] = 1.0;
    }
    else
    {
      dscale_[i] = 1.0 / ::sqrt( std::abs(dscale_[i]));
    }
  }

  matrix_->resetEvents ();

  events_ = 0;

  }


//-----------------------------------------------------------------------
//   connectToSolver_
//-----------------------------------------------------------------------


void AMGCLSolver::connectToSolver_ ()
{
}


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void AMGCLSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void AMGCLSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void AMGCLSolver::setEvents_ ( int events )
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


void AMGCLSolver::progressHandler_ ( idx_t jcol )
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


void AMGCLSolver::pivotHandler_ ( idx_t irow, double pivot )
{
}


//-----------------------------------------------------------------------
//   zeroPivotHandler_
//-----------------------------------------------------------------------


void AMGCLSolver::zeroPivotHandler_ ( idx_t irow, double pivot )
{
}


JIVE_END_PACKAGE( solver )

#endif