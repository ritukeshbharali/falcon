
/** @file cuDSSSolver.cpp
 *  @brief Wrapper to solve a linear system of equations with NVIDIA cuDSS.
 * 
 *  Requires: 
 *     - NVIDIA GPU with Compute Capability >=3.0
 *     - cuDSS installation
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 23 March 2024 
 *
 *  Updates (when, what and who)
 *     - [XX YYYYY 2024]
 *       
 */

#if defined(WITH_CUDSS)

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

#include "cuDSSSolver.h"


JEM_DEFINE_CLASS( jive::solver::cuDSSSolver );


JIVE_BEGIN_PACKAGE( solver )


using jem::System;
using jem::newInstance;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class cuDSSSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*   cuDSSSolver::TYPE_NAME          = "cuDSS";

const int     cuDSSSolver::NEW_VALUES_        = 1 << 0;
const int     cuDSSSolver::NEW_STRUCT_        = 1 << 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


cuDSSSolver::cuDSSSolver

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

  // cuDSS stuff

  cudaError_t cuda_err       = cudaSuccess;
  cudssStatus_t cudss_status = CUDSS_STATUS_SUCCESS;

  // Create cuda stream

  cuda_err = cudaStreamCreate(&stream);
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  // Create the cuDSS library handle

  cudss_status = cudssCreate(&handle);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssSetStream(handle, stream);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssConfigCreate(&solverConfig);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssDataCreate(handle, &solverData);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  mtype = CUDSS_MTYPE_GENERAL;
  mview = CUDSS_MVIEW_FULL;
  base  = CUDSS_BASE_ZERO;

  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


cuDSSSolver::~cuDSSSolver ()
{
  // Remember to destroy the GPU resources

  freeDeviceMemory_();

  cudssStatus_t cudss_status = CUDSS_STATUS_SUCCESS;

  cudss_status = cudssMatrixDestroy(A);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssMatrixDestroy(x);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssMatrixDestroy(b);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssDataDestroy(handle, solverData);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssConfigDestroy(solverConfig);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  cudss_status = cudssDestroy(handle);
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_(); 
}


//-----------------------------------------------------------------------
//   start
//-----------------------------------------------------------------------


void cuDSSSolver::start ()
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


void cuDSSSolver::finish ()
{
  if ( started_ )
  {
    started_--;
  }
}


//-----------------------------------------------------------------------
//   clear
//-----------------------------------------------------------------------


void cuDSSSolver::clear ()
{
  JEM_PRECHECK ( ! started_ );

  events_ = ~0x0;
}


//-----------------------------------------------------------------------
//   improve
//-----------------------------------------------------------------------


void cuDSSSolver::improve

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

  cudaError_t cuda_err       = cudaSuccess;
  cudssStatus_t cudss_status = CUDSS_STATUS_SUCCESS;

  while ( iiter_ < 1 )
  {
    // Upload lhs and rhs to GPU

    cuda_err = cudaMemcpy ( du_d, du.addr(), sizeof(double) * (dofCount), cudaMemcpyHostToDevice );
    if (cuda_err != cudaSuccess) freeDeviceMemory_();

    cuda_err = cudaMemcpy ( r_d, r.addr(), sizeof(double) * (dofCount), cudaMemcpyHostToDevice );
    if (cuda_err != cudaSuccess) freeDeviceMemory_();

    // Create the dense matrix on device for lhs and rhs

    cudss_status = cudssMatrixCreateDn(&b, 1, 1, dofCount, r_d, CUDA_R_64F,
                                       CUDSS_LAYOUT_COL_MAJOR);
    if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

    cudss_status = cudssMatrixCreateDn(&x, 1, 1, dofCount, du_d, CUDA_R_64F,
                                       CUDSS_LAYOUT_COL_MAJOR);
    if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

    // For some reason cuDSS factorization calls require the solution (du)
    // and the right hand side (r)

    // Symbolic factorization
    cudss_status = cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                                A, x, b);
    if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

    // Numeric factorization
    cudss_status = cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig, solverData,
                                A, x, b);
    if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

    // Solve
    cudss_status = cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
                                A, x, b);
    if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

    // Synchronize before returning control to host

    cuda_err = cudaStreamSynchronize(stream);
    if (cuda_err != cudaSuccess) freeDeviceMemory_();

    // Copy solution from device to host

    cuda_err = cudaMemcpy ( du.addr(), du_d, sizeof(double) * (dofCount), cudaMemcpyDeviceToHost );
    if (cuda_err != cudaSuccess) freeDeviceMemory_();

    // Update u += du

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


void cuDSSSolver::getInfo ( const Properties& info ) const
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


void cuDSSSolver::configure ( const Properties& props )
{}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void cuDSSSolver::getConfig ( const Properties& conf ) const
{}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void cuDSSSolver::setMode ( int mode )
{
  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int cuDSSSolver::getMode () const
{
  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void cuDSSSolver::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double cuDSSSolver::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* cuDSSSolver::getMatrix () const
{
  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* cuDSSSolver::getConstraints () const
{
  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   getNullSpace
//-----------------------------------------------------------------------


void cuDSSSolver::getNullSpace ( Matrix& nspace )
{
}


//-----------------------------------------------------------------------
//   setZeroThreshold
//-----------------------------------------------------------------------


void cuDSSSolver::setZeroThreshold ( double eps )
{
}


//-----------------------------------------------------------------------
//   getZeroThreshold
//-----------------------------------------------------------------------


double cuDSSSolver::getZeroThreshold () const
{
  return NAN;
}


//-----------------------------------------------------------------------
//   setMaxZeroPivots
//-----------------------------------------------------------------------


void cuDSSSolver::setMaxZeroPivots ( idx_t n )
{
}


//-----------------------------------------------------------------------
//   getMaxZeroPivots
//-----------------------------------------------------------------------


idx_t cuDSSSolver::getMaxZeroPivots () const
{
  return -1;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> cuDSSSolver::makeNew

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


void cuDSSSolver::declare ()
{
  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void cuDSSSolver::update_ ()
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

  // Allocate device memory for the matrix

  cudaError_t cuda_err       = cudaSuccess;
  cudssStatus_t cudss_status = CUDSS_STATUS_SUCCESS;

  cuda_err = cudaMalloc ( &offsets_d, sizeof(int)    * (msize + 1) );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  cuda_err = cudaMalloc ( &indices_d, sizeof(int)    * (nnz) );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  cuda_err = cudaMalloc ( &values_d,  sizeof(double) * (nnz) );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  // Upload matrix data to GPU

  cuda_err = cudaMemcpy ( offsets_d, offsets_.addr(), sizeof(int)    * (msize + 1), cudaMemcpyHostToDevice );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  cuda_err = cudaMemcpy ( indices_d, indices_.addr(), sizeof(int)    * (nnz), cudaMemcpyHostToDevice );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  cuda_err = cudaMemcpy ( values_d,  values_.addr(),  sizeof(double) * (nnz), cudaMemcpyHostToDevice );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  // Create the CSR matrix on device

  cudss_status = cudssMatrixCreateCsr( &A, msize, msize, nnz, offsets_d, NULL,
                                      indices_d, values_d, CUDA_R_32I, CUDA_R_64F,
                                      mtype, mview, base );
  if (cudss_status != CUDSS_STATUS_SUCCESS) freeDeviceMemory_();

  // Allocate device memory for the x and b

  cuda_err = cudaMalloc ( &du_d, sizeof(double) * (msize) );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  cuda_err = cudaMalloc ( &r_d, sizeof(double) * (msize) );
  if (cuda_err != cudaSuccess) freeDeviceMemory_();

  matrix_->resetEvents ();

  events_ = 0;

  }


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void cuDSSSolver::valuesChanged_ ()
{
  setEvents_ ( NEW_VALUES_ );
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void cuDSSSolver::structChanged_ ()
{
  setEvents_ ( NEW_STRUCT_ );
}


//-----------------------------------------------------------------------
//   setEvents_
//-----------------------------------------------------------------------


void cuDSSSolver::setEvents_ ( int events )
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
//   freeDeviceMemory_
//-----------------------------------------------------------------------


void cuDSSSolver::freeDeviceMemory_ ()
{
  cudaFree(offsets_d);
  cudaFree(indices_d);
  cudaFree(values_d);
  cudaFree(du_d);
  cudaFree(r_d);
}


JIVE_END_PACKAGE( solver )

#endif