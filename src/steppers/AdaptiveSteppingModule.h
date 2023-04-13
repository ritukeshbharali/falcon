
/** @file AdaptiveSteppingModule.h
 *  @brief Implements (time) stepping module with adaptive step-size.
 * 
 *  This class implements a module for Adaptive (time)
 *  stepping. If a certain choice of stepsize fails to
 *  converge, the stepsize is reduced. Again, if a step
 *  converges within 'optiter' iterations, the stepsize
 *  is increased. Once the stepsize is changed, the 
 *  module throws the action SolverNames::SET_STEP_SIZE.
 *  Models, such as the DirichletModel, modifies the
 *  stepsize.
 * 
 *  NOTE: This code is re-write of the Jive 2.2 version
 *        AdaptiveStepGeneralModule written by 
 *        F.P van der Meer,
 *        TU Delft, The Netherlands.
 *
 *  Date: May 2013
 * 
 *  Usage:      
 * 
 *    solver =
 *       {
 *        type       = "AdaptiveStepping";
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
 *  Updates (when, what and who)
 * 
 *  (2016) The original model had a default implementation for the Solver. The 
 *  modified implementation allows a solver to be given as input. This means 
 *  the method can be used for other solvers than the NonLin, even for implicit
 *  dynamics. (E.C. Simons, TU Delft) 
 *   
 *       
 */


#ifndef ADAPTIVE_STEP_GEN_MODULE_H
#define ADAPTIVE_STEP_GEN_MODULE_H

/* Include jem and jive headers */

#include <jem/io/PrintWriter.h>
#include <jem/util/CPUTimer.h>
#include <jive/app/Module.h>
#include <jive/model/Model.h>
#include <jive/implict/SolverModule.h>

using jem::Ref;
using jem::String;
using jem::NIL;
using jem::idx_t;
using jem::io::PrintWriter;
using jem::util::Properties;
using jem::util::CPUTimer;
using jive::app::Module;
using jive::model::Model;
using jive::implict::SolverModule;


//-----------------------------------------------------------------------
//   class AdaptiveSteppingModule
//-----------------------------------------------------------------------

/** @brief 
 *  The AdaptiveSteppingModule class implements an adaptive time-stepping
 *  method. 
 * 
 *  Acts as a solver module, embedded with a nonlinear solver.
 */

class AdaptiveSteppingModule : public Module
{
 public:

  typedef AdaptiveSteppingModule    Self;
  typedef Module                    Super;

  static const char*        TYPE_NAME;
  static const char*        SOLVER;
  static const char*        WRITE_STAT_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        START_INCR_PROP;
  static const char*        REDUCTION_PROP;
  static const char*        OPT_ITER_PROP;
  static const char*        STRICT_PROP;

  explicit                  AdaptiveSteppingModule

    ( const String&           name   = "AdaptiveStepping",
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

  virtual                  ~AdaptiveSteppingModule   ();

 private:

  bool                      reduceStep_

    ( const Properties&       globdat );

  bool                      increaseStep_

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

  double                    maxIncrTried_;

  bool                      writeStats_;

  Ref<PrintWriter>          statOut_;

  // for adaptive time stepping 

  double                    startIncr_;
  double                    minIncr_;
  double                    maxIncr_;
  double                    timeMax_;
  double                    reduction_;
  idx_t                     optIter_;

  double                    increment0_;
  double                    increment_;
  double                    timeA_;
  double                    timeB_;
  double                    time0_;
  double                    time_;

  // counters inside time step

  idx_t                     nCancels_;
  idx_t                     nContinues_;
  idx_t                     maxNIter_;
  bool                      triedSmall_;
  
  // new stuff by E.C. Simons
  bool                      strict_;

  // cumulative counters (for shutdown statistics)

  idx_t                     nRunTotal_;
  idx_t                     nCancTotal_;
  idx_t                     nContTotal_;
  idx_t                     nIterTotal_;
  idx_t                     nIterCont_;
  idx_t                     nCommTotal_;
  idx_t                     nCancCont_;      // canceled continues

  CPUTimer                  cpuTimer_;
  
};


#endif
