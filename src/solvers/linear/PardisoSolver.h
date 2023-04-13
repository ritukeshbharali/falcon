
/** @file PardisoSolver.h
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
#ifndef JIVE_SOLVER_PARDISO_H
#define JIVE_SOLVER_PARDISO_H

/* Include jem and jive headers */

#include <jem/base/Flags.h>
#include <jive/solver/import.h>
#include <jive/solver/DirectSolver.h>

#include <mkl_pardiso.h>
#include <mkl_service.h>


JIVE_BEGIN_PACKAGE( solver )


class Constrainer;


//-----------------------------------------------------------------------
//   Pardiso C prototype declaration
//-----------------------------------------------------------------------

extern "C" {

  void pardiso ( void* pt,  int* maxfct, int* mnum,  int* mtype,  int* phase,
                 int*  n,   double* a,   int* ia,    int* ja,     int* perm, 
                 int* nrhs, int* iparm,  int* msgvl, double* rhs, double* x,
                 int* error                                                  );
}

//-----------------------------------------------------------------------
//   class PardisoSolver
//-----------------------------------------------------------------------

/** @brief 
 *  The PardisoSolver class implements a wrapper for using the Pardiso
 *  direct solver from the Intel MKL library. 
 * 
 *  @note Does not work with MPI.
 */ 

class PardisoSolver : public DirectSolver
{
 public:

  JEM_DECLARE_CLASS       ( PardisoSolver, DirectSolver );

  static const char*        TYPE_NAME;
  static const double       PIVOT_THRESHOLD;
  static const int          MAX_ITER;
  static const char*        REORDER_METHODS[3];

  static const char*        MATRIX_TYPE;
  static const char*        PRINT_INFO;
  static const char*        NUM_THREADS;
  static const char*        SORT_COLUMNS;
  static const char*        MATRIX_CHECKER;  

  enum                      Option
  {
                              PRINT_PIVOTS = 1 << 0
  };

  typedef
    jem::Flags<Option>      Options;

  enum                      ReorderMethod
  {
                              REORDER_NONE,
                              REORDER_MATRIX,
                              REORDER_COLUMNS
  };


                            PardisoSolver

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

  void                      setPivotThreshold

    ( double                  alpha );

  double                    getPivotThreshold () const;
  Options                   getOptions        () const;

  void                      setOption

    ( Option                  option,
      bool                    yesno = true );

  void                      setOptions

    ( Options                 options );

  ReorderMethod             getReorderMethod  () const;

  void                      setReorderMethod

    ( ReorderMethod           method );

  static Ref<Solver>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       params,
      const Properties&       globdat );

  static void               declare           ();


 protected:

  virtual                  ~PardisoSolver          ();


 private:

  void                      update_             ();
  void                      connectToSolver_    ();
  void                      valuesChanged_      ();
  void                      structChanged_      ();

  void                      freePardisoMemory_  ();

  void                      setEvents_

    ( int                     events );

  void                      progressHandler_

    ( idx_t                   jcol  );

  void                      pivotHandler_

    ( idx_t                   irow,
      double                  pivot );

  void                      zeroPivotHandler_

    ( idx_t                   irow,
      double                  pivot );


 private:

  static const int          NEW_VALUES_;
  static const int          NEW_STRUCT_;


  Ref<AbstractMatrix>       matrix_;
  Ref<Constrainer>          conman_;
  Ref<Writer>               debug_;

  int                       mode_;
  double                    small_;
  double                    smallPiv_;
  double                    precision_;
  ReorderMethod             reorder_;
  Options                   options_;
  idx_t                     maxZeroes_;
  bool                      updated_;

  int                       iiter_;
  double                    error_;
  int                       events_;
  idx_t                     started_;

  void*                     pt_[64];
  int                       iparm_[64];

  int                       n_;
  int                       nrhs_;
  int                       maxfct_;
  int                       mnum_;
  int                       mtype_;
  int                       phase_;
  int                       msglvl_;
  int                       sortCols_;
  int                       matrixChk_;

  int                       idum_;
  double                    ddum_;
  int                       errorFlag_;
  int                       nThreads_;

  jem::Array<int>           offsets_;
  jem::Array<int>           indices_;
  Vector                    values_;
  

};


JEM_DEFINE_FLAG_OPS( PardisoSolver::Options )


JIVE_END_PACKAGE( solver )

#endif
