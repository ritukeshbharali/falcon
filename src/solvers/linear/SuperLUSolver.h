
/** @file SuperLUSolver.h
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

#ifndef JIVE_SOLVER_AMGCL_H
#define JIVE_SOLVER_AMGCL_H

/* Include jem and jive headers */

#include <jem/base/Flags.h>
#include <jive/solver/import.h>
#include <jive/solver/DirectSolver.h>

#include <slu_mt_ddefs.h>
#include <vector>

JIVE_BEGIN_PACKAGE( solver )


class Constrainer;

//-----------------------------------------------------------------------
//   class SuperLUSolver
//-----------------------------------------------------------------------

/** @brief 
 *  The \c SuperLUSolver class implements a wrapper for solving
 *  the system of linear equations using the multi-thread version
 *  of the SuperLU library (SuperLU_MT). SuperLU operates on 
 *  Compressed Sparse Column (CSC) matrix opposed as the Jive 
 *  Compressed Sparse Row (CSR) matrix. Therefore, we solve
 *  A^T X = B.
 * 
 *  \note Does not work with MPI.
 * 
 *  Below is an example how the linear SuperLU solver is defined
 *  with a nonlinear solver:
 *  
 *  \code
 *  solver = "Nonlin"
    {
      precision  = 1.e-4;
      maxIter    = 100;
      lineSearch = true;
      solver =
      {
        type = "SuperLU";   // SuperLU
        numThreads = 2;     // Number of threads
      };
    };
 *  \endcode
 *   
 */ 

class SuperLUSolver : public DirectSolver
{
 public:

  JEM_DECLARE_CLASS       ( SuperLUSolver, DirectSolver );

  static const char*        TYPE_NAME;
  static const double       PIVOT_THRESHOLD;
  static const int          MAX_ITER;
  static const char*        REORDER_METHODS[3];

  static const char*        NUM_THREADS;

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


                            SuperLUSolver

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

  virtual                  ~SuperLUSolver          ();


 private:

  void                      update_           ();
  void                      connectToSolver_  ();
  void                      valuesChanged_    ();
  void                      structChanged_    ();


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

  int                       iiter_;
  double                    error_;
  int                       events_;
  idx_t                     started_;

  // CSR arrays

  jem::Array<int>           offsets_;    //!< CSR row offsets
  jem::Array<int>           jindices_;   //!< CSR col indices
  Vector                    values_;     //!< CSR values

  // SuperLU parameters

  int                       nThreads_;   //!< number of threads
  int                       cPermSpec_;  //!< col. perm. type
  int                       panelSize_;  //!< cols treated as unit task
  int                       relax_;      //!< cols grouped as relaxed supernode
  yes_no_t                  refact_;     //!< re-factor
  yes_no_t                  usePr_;
  fact_t                    fact_;       //!< factorized or not
  trans_t                   trans_;      //!< solve A^T X = B
  double                    diagPivThr_; //!< diagonal pivot threshold
  double                    dropTol_;    //!< drop tolerance
  void*                     work_;       //!< pre-allocated work space
  int                       lwork_;      //!< length of work

  superlumt_options_t       sluOpts_;    //!< solver options
  Gstat_t                   Gstat_;      //!< solver stats
  bool                      updated_;    //!< update flag
  int                       info_;       //!< info

  // SuperLU matrices and vectors

  SuperMatrix               A_;          //!< Matrix
  SuperMatrix               AC_;         //!< Permuted matrix A
  SuperMatrix               L_;          //!< L matrix of A
  SuperMatrix               U_;          //!< U matrix of A
  SuperMatrix               B_;          //!< rhs 
  SuperMatrix               X_;          //!< lhs (solution)
  std::vector<int>          cPerm_;      //!< col. perm
  std::vector<int>          rPerm_;      //!< row. perm

};


JEM_DEFINE_FLAG_OPS( SuperLUSolver::Options )


JIVE_END_PACKAGE( solver )

#endif
#endif
