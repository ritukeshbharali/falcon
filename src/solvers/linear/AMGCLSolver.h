
/** @file AMGCLSolver.h
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

#ifndef JIVE_SOLVER_AMGCL_H
#define JIVE_SOLVER_AMGCL_H

/* Include jem and jive headers */

#include <jem/base/Flags.h>
#include <jive/solver/import.h>
#include <jive/solver/DirectSolver.h>

#include<vector>

JIVE_BEGIN_PACKAGE( solver )


class Constrainer;

//-----------------------------------------------------------------------
//   class AMGCLSolver
//-----------------------------------------------------------------------

/** @brief 
 *  The AMGCLSolver class implements a wrapper for using the ANGCL
 *  library solvers and preconditioners.
 * 
 *  @note Does not work with MPI.
 */ 

class AMGCLSolver : public DirectSolver
{
 public:

  JEM_DECLARE_CLASS       ( AMGCLSolver, DirectSolver );

  static const char*        TYPE_NAME;
  static const double       PIVOT_THRESHOLD;
  static const int          MAX_ITER;
  static const char*        REORDER_METHODS[3];

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


                            AMGCLSolver

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

  virtual                  ~AMGCLSolver          ();


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

  jem::Array<int>           offsets_;    // row offsets
  jem::Array<int>           jindices_;   // col indices
  Vector                    values_;     // values
  Vector                    dscale_;     // diagonal entries

};


JEM_DEFINE_FLAG_OPS( AMGCLSolver::Options )


JIVE_END_PACKAGE( solver )

#endif
#endif
