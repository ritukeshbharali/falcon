
/** @file AMGXSolver.h
 *  @brief Wrapper to solve a linear system with AMGX.
 * 
 *  Requires: 
 *     - NVIDIA GPU with Compute Capability >=3.0
 *     - AMGX installation with dependencies
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 02 January 2024 
 *
 *  Updates (when, what and who)
 *     - [XX YYYYY 2024]
 *       
 */

#if defined(WITH_AMGX)

#ifndef JIVE_SOLVER_AMGX_H
#define JIVE_SOLVER_AMGX_H

/* Include jive headers */

#include <jive/solver/import.h>
#include <jive/solver/DirectSolver.h>

/* Include AMGX header file */

#include "amgx_c.h"

JIVE_BEGIN_PACKAGE( solver )


class Constrainer;


//=======================================================================
//   class AMGXSolver
//=======================================================================

/** @brief Implements NVIDIA's Algebraic multi-grid (AmgX) Solver
 * 
 *  The class \c AMGXSolver implements a wrapper to solve the linear 
 *  system of equations using the Algebraic Multi-grid (AmgX) solver,
 *  distributed by NVIDIA <a href="https://github.com/NVIDIA/AMGX" target="_blank">(Link to Github repo)</a>.
 * 
 *  \note Wrapper works only for sequential (shared-memory) runs. 
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
        type = "AmgX"; 
      };
    };
 *  \endcode
 * 
 *  Additionally, a 'solver.json' file is required that configures
 *  the AmgX solver options. See <a href="https://github.com/NVIDIA/AMGX/tree/main/src/configs" target="_blank">(Github repo directory)</a>
 *  for inspiration.  
 * 
 */ 

class AMGXSolver : public DirectSolver
{
 public:

  JEM_DECLARE_CLASS       ( AMGXSolver, DirectSolver );

  static const char*        TYPE_NAME;

                            AMGXSolver

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

  virtual                  ~AMGXSolver          ();


 private:

  void                      update_           ();
  void                      valuesChanged_    ();
  void                      structChanged_    ();

  void                      setEvents_

    ( int                     events );


 private:

  static const int          NEW_VALUES_;
  static const int          NEW_STRUCT_;


  Ref<AbstractMatrix>       matrix_;
  Ref<Constrainer>          conman_;
  Ref<Writer>               debug_;

  int                       mode_;
  double                    precision_;

  int                       iiter_;
  double                    error_;
  int                       events_;
  idx_t                     started_;

  jem::Array<int>           offsets_;
  jem::Array<int>           indices_;
  Vector                    values_;

  // AMGX Data structures

  AMGX_config_handle        AMGXcfg_;
  AMGX_resources_handle     AMGXrsrc_;
  AMGX_matrix_handle        AMGXA_;
  AMGX_vector_handle        AMGXb_, AMGXx_;
  AMGX_solver_handle        AMGXsolver_;
  AMGX_SOLVE_STATUS         AMGXstatus_;

};


JIVE_END_PACKAGE( solver )

#endif
#endif