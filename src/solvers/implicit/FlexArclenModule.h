/** 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implements a flexible path following
 *  method using both displacement control and energy based 
 *  arclength control. This is very helpful in tracing
 *  complex equilibrium paths.
 *
 *  Basic idea:
 *
 *    Start the simulation with load control until dissipation 
 *    becomes significant or convergence becomes problematic.
 *    Then continue with energy-based arclength
 *    control with amount of released energy computed. Proceeding
 *    with this arclength control till hardening branch is 
 *    detected, then switch back to displacement control.
 *  
 *  Authors: V.P. Nguyen, F.P. van der Meer
 *  Date: January 2009
 *
 *  Updates (when, what and who)
 *     - [10 January 2024] Added the time-step computing
 *       arclength module. (RB)
 * 
 */

#ifndef FLEX_ARCLEN_MODULE_H
#define FLEX_ARCLEN_MODULE_H

#include <jem/io/PrintWriter.h>
#include <jem/util/CPUTimer.h>
#include <jive/app/Module.h>
#include <jive/model/Model.h>

namespace jive
{
  namespace implict
  {
    class ArclenModule;
    class TSArclenModule;
    class NonlinModule;
  }
}

using jem::Ref;
using jem::String;
using jem::NIL;
using jem::idx_t;
using jem::lint;
using jem::io::PrintWriter;
using jem::util::Properties;
using jem::util::CPUTimer;
using jive::app::Module;
using jive::model::Model;
using jive::implict::NonlinModule;
using jive::implict::ArclenModule;
using jive::implict::TSArclenModule;
using jive::implict::SolverModule;


//-----------------------------------------------------------------------
//   class FlexArclenModule
//-----------------------------------------------------------------------


class FlexArclenModule : public Module
{
 public:

  typedef FlexArclenModule   Self;
  typedef Module             Super;

  static const char*        TYPE_NAME;
  static const char*        TIMESTEP_BASED;
  static const char*        NONLIN;
  static const char*        ARCLEN;
  static const char*        TSARCLEN;
  static const char*        TMP_NONLIN_PROP;
  static const char*        KEEP_LARGE_PROP;
  static const char*        WRITE_STAT_PROP;

  explicit                  FlexArclenModule

    ( const String&         name    = "FlexArclen",
      Ref<NonlinModule>     solver1 = NIL,
      Ref<ArclenModule>     solver2 = NIL,
      Ref<TSArclenModule>   solver3 = NIL );

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
      const Properties&       globdat )    const;

  static Ref<Module>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

 protected:

  virtual                  ~FlexArclenModule  ();

 private:

  bool                      switchSolver_

    ( const Properties&       globdat );

  bool                      switchToNonlin_

    ( const Properties&       globdat,
      const Properties&       params = Properties() );

  bool                      switchToArclen_

    ( const Properties&       globdat,
      const Properties&       params = Properties() );

  bool                      reduceStep_

    ( const Properties&       globdat );

  bool                      increaseStep_

    ( const Properties&       globdat );

  bool                      turnOffResSwt_

    ( const Properties&       globdat );

  bool                      considerCarefulMode_

    ( const Properties&       globdat );

  void                      storeMaxIncrs_

    ( const Properties&       globdat );

  void                      commit_ 
    
    ( const Properties&       globdat );

  void                      continue_ 
    
    ( const Properties&       params,
      const Properties&       globdat );

  void                      cancel_
    
    ( const Properties&       globdat );

  void                      getDissForce_

    ( const Properties&       globdat ) const;

 private:

  Ref<NonlinModule>         nonlinSolver_;     // std. Newton Raphson
  Ref<ArclenModule>         arclenSolver_;     // load-disp Arclen
  Ref<TSArclenModule>       tsArclenSolver_;   // time-step Arclen

  Ref<SolverModule>         currentSolver_;

  Ref<Model>                model_;

  idx_t                     istep_;
  idx_t                     istep0_;

  bool                      tsBased_;

  double                    maxDispTried_;
  double                    maxArcTried_;
  double                    swtResidual0_;
  double                    arcBackup_;

  bool                      nonlinTriedAll_;
  bool                      arclenTriedAll_;
  bool                      arclenTriedSma_;
  bool                      persistent_;     // don't switch until cancel/commit
  bool                      writeStats_;

  Ref<PrintWriter>          statOut_;

  // counters inside time step

  idx_t                     nCancels_;
  idx_t                     nContinues_;

  // cumulative counters (for shutdown statistics)

  lint                      nRunTotal_;
  lint                      nNonlTotal_;
  lint                      nCancTotal_;
  lint                      nContTotal_;
  lint                      nIterTotal_;
  lint                      nIterCont_;
  lint                      nCommTotal_;
  lint                      nCancCont_;      // canceled continues
  lint                      nNonlCanc_;      // cancels in nonlin
  lint                      nChanges_;

  CPUTimer                  cpuTimer_;

  // some flags for trying to switch to nonlin after crack propagation

  bool                      doTmpNonlin_;    // general flag from input
  bool                      isTmpNonlin_;    // current state
  bool                      keepLargeTmp_;   // keep solution when large
                                             // dissipation incr in tmpNonlin
  idx_t                     doneTmpNonlin_;  // do only once (per step size)
                                             // 0: not done
                                             // 1: doing
                                             // 2: done and turned off
  bool                      triedCareful_;
  bool                      isCareful_;

  // some flags to temporarily turn of residual based switch in XArlenModule

  bool                      noResSwitch_;
  bool                      doneResSwitch_;
};


#endif
