
/** @file ReduceSteppingModule.h
 *  @brief Implements (time) stepping module with decreasing step-size.
 * 
 *  This class implements a module for reducing (time) 
 *  step sizes based on user input or when the solver
 *  fails to converge. Unlike the Adaptive Stepping 
 *  Module, the step size is never increased based on
 *  optimal iterations per step. Once the stepsize is 
 *  changed, the module throws the action 
 *  SolverNames::SET_STEP_SIZE. Models, such as the 
 *  DirichletModel, modifies the stepsize.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: XX March 2022
 * 
 *  NOTE: This code is based on the Adaptive Step
 *        GeneralModule written by F.P van der Meer,
 *        TU Delft, The Netherlands.
 * 
 *  Usage: ( typically for brittle failure )
 * 
 *    solver =
 *       {
 *        type   = "ReduceStepping";
 *        
 *        nonlin = 
 *           {
 *             precision   = 1.e-6;
 *             solver.type = "SparseLU";
 *             maxIter     = 10; 
 *           }
 *        optIter   = 4;
 *        reduction = 0.5;
 *        startIncr = 1.e+0;
 *        minIncr   = 1.e-5;
 *        maxIncr   = 1.e+2;
 *        strict    = true;
 *
 *       };
 *       
 */


#ifndef REDUCE_STEPPING_MODULE_H
#define REDUCE_STEPPING_MODULE_H

/* Include jem and jive headers */

#include <jive/app/Module.h>
#include <jive/model/Model.h>
#include <jive/implict/SolverModule.h>

using jem::Ref;
using jem::String;
using jem::NIL;
using jem::idx_t;
using jem::util::Properties;
using jem::util::CPUTimer;
using jive::app::Module;
using jive::model::Model;
using jive::implict::SolverModule;


//-----------------------------------------------------------------------
//   class ReduceSteppingModule
//-----------------------------------------------------------------------

/** @brief 
 *  The ReduceSteppingModule class implements an adaptively decreasing 
 *  time-stepping method. 
 * 
 *  Acts as a solver module, embedded with a nonlinear solver. Used
 *  mainly for brittle fracture models.
 */

class ReduceSteppingModule : public Module
{
 public:

  typedef ReduceSteppingModule    Self;
  typedef Module                  Super;

  static const char*        TYPE_NAME;
  static const char*        SOLVER;
  static const char*        START_INCR_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        REDUCE_STEP_PROP;
  static const char*        REDUCTION_PROP;
  
  explicit                  ReduceSteppingModule

    ( const String&           name   = "ReduceStepping",
      Ref<SolverModule>       solver = NIL                 );

  virtual Status            init

    ( const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual Status            run

    ( const Properties&       globdat );

  virtual void              shutdown

    ( const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       props,
      const Properties&       globdat )        const;

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       globdat,
      const Properties&       props );

 protected:

  virtual                  ~ReduceSteppingModule   ();

 private:

  bool                      reduceStep_

    ( const Properties&       globdat );

  void                      commit_ 
    
    ( const Properties&       globdat );

  void                      continue_ 
    
    ( const Properties&       params,
      const Properties&       globdat );

  void                      cancel_
    
    ( const Properties&       globdat );

  void                      setStepSize_

    ( const Properties&       globdat ) const;

 private:

  Ref<SolverModule>         solver_;

  Ref<Model>                model_;

  idx_t                     istep_;
  idx_t                     istep0_;

  // for time stepping 

  double                    startIncr_;
  double                    minIncr_;
  int                       cutStep_;
  double                    reduction_;
  
  // stores the current and old time step size

  double                    increment0_;
  double                    increment_;
  
};


#endif
