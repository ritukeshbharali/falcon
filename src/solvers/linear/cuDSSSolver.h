
/** @file cuDSSSolver.h
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

#ifndef JIVE_SOLVER_CUDSS_H
#define JIVE_SOLVER_CUDSS_H

/* Include jive headers */

#include <jive/solver/import.h>
#include <jive/solver/DirectSolver.h>

/* Include cuda headers */

#include <cuda_runtime.h>
#include "cuDSS.h"

JIVE_BEGIN_PACKAGE( solver )  

class Constrainer;


//=======================================================================
//   class cuDSSSolver
//=======================================================================

/** @brief Implements NVIDIA's GPU-accelerated Direct Sparse Solver
 * 
 *  The class \c cuDSSSolver implements a wrapper to solve the linear 
 *  system of equations using the cuDSS solver, which is NVIDIA's
 *  optimized, first-generation GPU-accelerated Direct Sparse Solver 
 *  library.
 * 
 *  \note Works only for sequential (shared-memory) runs. 
 * 
 *  \todo Add configuration options.
 * 
 *  Below is an usage example, demonstrating how the solver is defined in
 *  the input file:
 * 
 *  \code
 *  solver = 
    {
      type = "Nonlin";
    
      precision = 1.e-6;
    
      maxIter   = 5;
  
      solver =
      {
        type = "cuDSS"; 
      };
    };
 *  \endcode
 * 
 */ 

class cuDSSSolver : public DirectSolver
{
 public:

  JEM_DECLARE_CLASS       ( cuDSSSolver, DirectSolver );

  static const char*        TYPE_NAME;

                            cuDSSSolver

    ( const String&           name,
      Ref<AbstractMatrix>     matrix,
      Ref<Constraints>        cons );

  virtual void              start             ()       override;
  virtual void              finish            ()       override;
  virtual void              clear             ()       override;

  virtual void              improve

    ( const Vector&           lhs,
      const Vector&           rhs )                    override;

  virtual void              getInfo

    ( const Properties&       info )             const override;

  virtual void              configure

    ( const Properties&       props )                  override;

  virtual void              getConfig

    ( const Properties&       props )            const override;

  virtual void              setMode

    ( int                     mode )                   override;

  virtual int               getMode           () const override;

  virtual void              setPrecision

    ( double                  eps )                    override;

  virtual double            getPrecision      () const override;
  virtual AbstractMatrix*   getMatrix         () const override;
  virtual Constraints*      getConstraints    () const override;

  virtual void              getNullSpace

    ( Matrix&                 nspace )                 override;

  virtual void              setZeroThreshold

    ( double                  eps )                    override;

  virtual double            getZeroThreshold  () const override;

  virtual void              setMaxZeroPivots

    ( idx_t                   maxPivots )              override;

  virtual idx_t             getMaxZeroPivots  () const override;

  static Ref<Solver>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       params,
      const Properties&       globdat );

  static void               declare           ();


 protected:

  virtual                  ~cuDSSSolver          ();


 private:

  void                      update_           ();
  void                      valuesChanged_    ();
  void                      structChanged_    ();

  void                      freeDeviceMemory_ ();

  void                      setEvents_

    ( int                     events );


 private:

  static const int          NEW_VALUES_;
  static const int          NEW_STRUCT_;


  Ref<AbstractMatrix>       matrix_;       //!< Pointer to the Jive matrix
  Ref<Constrainer>          conman_;       //!< Pointer to the Jive constrainer
  Ref<Writer>               debug_;        //!< Pointer to the Jive printer

  int                       mode_;
  double                    precision_;    //!< Precision

  int                       iiter_;        //!< Iteration counter
  double                    error_;        //!< Error norm(Ax-b)
  int                       events_;
  idx_t                     started_;

  jem::Array<int>           offsets_;      //!< Row offsets of CSR matrix
  jem::Array<int>           indices_;      //!< Column indices of CSR matrix
  Vector                    values_;       //!< Values of CSR matrix

  // cuDSS Data structures

  cudssMatrix_t A, x, b;                   //!< cudss Matrices for A,x and b
                                           //!< A is of type CSR, x and b are Dense

  cudssMatrixType_t     mtype;  // CUDSS_MTYPE_GENERAL
  cudssMatrixViewType_t mview;  // CUDSS_MVIEW_FULL
  cudssIndexBase_t      base;   // CUDSS_BASE_ZERO

  cudaStream_t          stream = NULL;
  cudssHandle_t         handle;
  cudssConfig_t         solverConfig;
  cudssData_t           solverData;

  int    *offsets_d;    //!< Pointer to row offsets on device
  int    *indices_d;    //!< Pointer to column indices on device
  double *values_d;     //!< Pointer to values on device
  double *du_d;         //!< Pointer to solution (x) on device
  double *r_d;          //!< Pointer to residual (Ax-b) on device

};


JIVE_END_PACKAGE( solver )

#endif
#endif