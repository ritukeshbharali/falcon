
/** @file UmfpackSolver.h
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


#ifndef JIVE_SOLVER_UMFPACK_H
#define JIVE_SOLVER_UMFPACK_H

/* Include jem and jive headers */

#include <jem/base/Flags.h>
#include <jive/solver/import.h>
#include <jive/solver/DirectSolver.h>
#include <umfpack.h>
#include <SuiteSparse_config.h>


JIVE_BEGIN_PACKAGE( solver )


class Constrainer;


//-----------------------------------------------------------------------
//   class UmfpackSolver
//-----------------------------------------------------------------------

/** @brief 
 *  The UmfpackSolver class implements a wrapper for using the Umfpack
 *  direct solver from SuiteSparse. 
 * 
 *  @note Does not work with MPI.
 */ 

class UmfpackSolver : public DirectSolver
{
 public:

  JEM_DECLARE_CLASS       ( UmfpackSolver, DirectSolver );

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


                            UmfpackSolver

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

  virtual                  ~UmfpackSolver          ();


 private:

  void                      update_           ();
  void                      connectToSolver_  ();
  void                      valuesChanged_    ();
  void                      structChanged_    ();

  void                      freeSymbolic_     ();
  void                      freeNumeric_      ();

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

  void*                     symbolic_;
  void*                     numeric_;
  double                    control_[UMFPACK_CONTROL];
  double                    info_[UMFPACK_INFO];
  jem::Array<long>          offsets_;
  jem::Array<long>          indices_;
  Vector                    values_;
  bool                      updated_;

};


JEM_DEFINE_FLAG_OPS( UmfpackSolver::Options )


JIVE_END_PACKAGE( solver )

#endif
