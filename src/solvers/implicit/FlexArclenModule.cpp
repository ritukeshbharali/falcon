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

#include <jem/base/limits.h>
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/base/Exception.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/io/FileWriter.h>
#include <jem/io/NumberFormat.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/util/Globdat.h>
#include <jive/util/utilities.h>
#include <jive/util/DofSpace.h>
#include <jive/implict/SolverInfo.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/ArclenModule.h>
#include <jive/app/ModuleFactory.h>

#include "TSArclenModule.h"
#include "FlexArclenModule.h"
#include "util/XNames.h"


using jem::IllegalOperationException;
using jem::System;
using jem::max;
using jem::maxOf;
using jem::newInstance;
using jem::Error;
using jem::io::endl;
using jem::io::FileWriter;
using jem::io::NumberFormat;
using jem::util::StringUtils;
using jive::Vector;
using jive::implict::SolverInfo;
using jive::util::Globdat;
using jive::util::joinNames;
using jive::util::DofSpace;



//=======================================================================
//   class FlexArclenModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* FlexArclenModule::TYPE_NAME              = "FlexArclen";
const char* FlexArclenModule::TIMESTEP_BASED         = "timeStepBased";
const char* FlexArclenModule::TMP_NONLIN_PROP        = "doTmpNonlin";
const char* FlexArclenModule::NONLIN                 = "nonLin";
const char* FlexArclenModule::ARCLEN                 = "arcLen";
const char* FlexArclenModule::TSARCLEN               = "tsArcLen";
const char* FlexArclenModule::KEEP_LARGE_PROP        = "keepLargeTmp";
const char* FlexArclenModule::WRITE_STAT_PROP        = "writeStats";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

FlexArclenModule::FlexArclenModule 

  ( const String&        name,
    Ref<NonlinModule>    solver1,
    Ref<ArclenModule>    solver2,
    Ref<TSArclenModule>  solver3 ) :

      Super           ( name    ),
      nonlinSolver_   ( solver1 ),
      arclenSolver_   ( solver2 ),
      tsArclenSolver_ ( solver3 ),
      currentSolver_  ( solver1 ),

      istep_          ( 0 ),
      istep0_         ( 0 ),
      nCancels_       ( 0 ),
      nContinues_     ( 0 ),
      nRunTotal_      ( 0 ),
      nCancTotal_     ( 0 ),
      nContTotal_     ( 0 ),
      nIterTotal_     ( 0 ),
      nCommTotal_     ( 0 ),
      nCancCont_      ( 0 ),
      nChanges_       ( 0 ),
      doneTmpNonlin_  ( 0 )


{
  maxDispTried_   = 0.;
  maxArcTried_    = 0.;
  swtResidual0_   = 0.;
  arcBackup_      = 0.;

  writeStats_     = false;
  doTmpNonlin_    = false;
  keepLargeTmp_   = true;
  isTmpNonlin_    = false;
  nonlinTriedAll_ = false;
  arclenTriedAll_ = false;
  arclenTriedSma_ = false;
  persistent_     = false;
  noResSwitch_    = false;
  doneResSwitch_  = false;
  triedCareful_   = false;

  tsBased_        = false;
}

FlexArclenModule::~FlexArclenModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status FlexArclenModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  model_ = Model::get ( globdat, getContext() );
  istep_ = istep0_ = 0;

  globdat.set ( Globdat::TIME_STEP, istep_ );
  globdat.set ( "var.accepted", false );

  nonlinSolver_  ->init ( conf, props, globdat );
  arclenSolver_  ->init ( conf, props, globdat );
  tsArclenSolver_->init ( conf, props, globdat );

  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status FlexArclenModule::run ( const Properties& globdat )

{
  using jem::Exception;

  Properties  info      = SolverInfo::get ( globdat );
  Properties  params;

  bool        accept    = true;
  bool        discard   = false;
  bool        converged = false;
  String      terminate;
  String      doSwitch;

  nRunTotal_  ++;
  nNonlTotal_ += ( currentSolver_ == nonlinSolver_ );

  info.clear ();

  currentSolver_->advance ( globdat );

  // Store maximum values of increment tried in this time step

  if ( nContinues_ == 0 )
  {
    // storeMaxIncrs_ ( globdat );
  }
  
  // Advance the time/load step number

  istep_ = istep0_;

  if ( istep_ < maxOf( istep_ ) )
  {
    istep_++;
  }

  globdat.set ( Globdat::TIME_STEP, istep_ );

  // try to find a solution

  try
  {
    currentSolver_->solve ( info, globdat );

    converged = true;
  }
  catch ( const jem::Exception& ex )
  {
    System::out() << "FlexArclenModule caught a solver exception: \n" 
                  << "   " << ex.where() << "\n"
                  << "   " << ex.what() << "\n\n";

    discard    = true;
  }

  // get # of iterations, NB: not set when error thrown by solver
  
  idx_t          iterCount = 0; 
  info.find   (  iterCount , SolverInfo::ITER_COUNT );

  nIterTotal_ += iterCount;

  if ( nContinues_ > 0 ) nIterCont_ += iterCount;

  if ( converged )
  {
    // Ask the model whether to accept converged solution
    // (Model can ask for switch via DO_SWITCH and for cancel via DISCARD)

    idx_t newCh = 0;

    params.set( SolverInfo::ITER_COUNT, iterCount       );
    params.set( XActions::CHECK_BOUNDS, ! keepLargeTmp_ );

    model_->takeAction ( XActions::CHECK_COMMIT, params, globdat );

    params.find ( discard, XActions::DISCARD      );
    params.find ( accept,  XActions::ACCEPT       );
    params.find ( newCh,   XActions::CHANGE_COUNT );

    nChanges_ += newCh;
  }

  // set flag in globdat for output indicating solution quality

  globdat.set ( "var.accepted", ( accept && ! discard ) );

  if ( ! discard )
  {
    // in CHECK_COMMIT, a check for switch has also been done:
    // do switch if the model has asked for it

    params.find ( doSwitch, XActions::DO_SWITCH );

    if ( ( doSwitch == "please"   && ! persistent_ ) ||
         ( doSwitch == "ifAccept" && accept        ) )
    {
      switchSolver_ ( globdat );
    }

    // take measures based on result of CHECK_COMMIT

    if ( accept )
    {
      // The final solution for this time step has been obtained
      // (possibly switch from tmpNonlin)

      commit_ ( globdat );

      if ( params.find ( terminate, XActions::TERMINATE ) ) 
      {
        if ( terminate == "sure" )
        {
          System::out() << "Some model gave the terminate signal during "
            << "checkCommit,\nFlexArclenModule terminates run.\n"
            << "NB: output modules are not called for this time step\n";
          return EXIT;
        }
      }
    }
    else
    {
      // We have an equilibrium solution but it is not accepted:
      // continue with this time step
      // (possibly switch to tmpNonlin)

      continue_ ( params, globdat );
    }
  }
  else
  {
    // 1. discard solution

    cancel_ ( globdat );

    // 2. find proper strategy to try again
   
    // ****
    // if the TimeStepArclenModule has quit because high residual scale factor
    // --> switch to nonlin

    if ( iterCount == 1 && nContinues_ == 0 && switchToNonlin_ ( globdat ) ) 
    {
      System::out() << "resSwitch: Returning to the same load step with "
                    << currentSolver_->getName() << "\n";

      doneResSwitch_  = true;

      // supposedly smaller energy increments will not work either

      arclenTriedAll_ = true;
    }

    // ****
    // switch to careful mode

    if ( nContinues_ == 1 && nChanges_ == 1 && !triedCareful_ &&
         considerCarefulMode_(globdat) )
    {
      System::out() << "Returning to the same load step in careful mode\n";
    }
    
    // ****
    // if tmpNonlin was tried, try again without tmpNonlin 
    // (same step size!)
       
    else if ( doneTmpNonlin_ == 1 )
    {
      ++doneTmpNonlin_;

      isTmpNonlin_   = false;

      switchToArclen_ ( globdat );

      Properties backup;

      backup.set ( XProps::STEP_SIZE, arcBackup_ );

      model_->takeAction ( XActions::SET_STEP_SIZE, backup, globdat );

      System::out() << "Returning to the same load step "
                    << "with tmpNonlin option turned off\n";
    }

    // ****
    // try to reduce the step size

    else if ( reduceStep_ ( globdat ) )
    {
      System::out() << "Returning to the same load step "
                    << "with reduced step\n";
      doneTmpNonlin_ = 0;
    }

    // ****
    // switch solver

    else if ( switchSolver_ ( globdat ) )
    {
      System::out() << "Returning to the same load step, now with "
                    << currentSolver_->getName() << "\n";
    }

    // ****
    // try arclength control with increasing step size

    else if ( increaseStep_ ( globdat ) )
    {
      System::out() << "Returning to the same load step "
                    << "with increased step\n";
    }

    // ****
    // try again without residual based switch

    else if ( turnOffResSwt_ ( globdat ) )
    {
      System::out() << "Returning tot the same load step,\n " 
                    << "now without residual based switch in TimeStepArclenModule\n";
    }

    // ****
    // give it up

    else
    {
      System::out() << "\n*** FlexArclenModule is out of inspiration,\n"
                    << "    can't find a winning strategy for this step."
                    << " Sorry.\n\n";

      return EXIT;
    }
    nContinues_ = 0;
    nChanges_   = 0;
  }

  // Always reassemble dissipation force. Is not always necessary, 
  // especially after cancel without continues, but ThermalModule might 
  // add it's part anyway

  getDissForce_ ( globdat );  

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void FlexArclenModule::shutdown ( const Properties& globdat )
{
  nonlinSolver_  ->shutdown ( globdat );
  arclenSolver_  ->shutdown ( globdat );
  tsArclenSolver_->shutdown ( globdat );

  System::out() << "FlexArclenModule statistics ..." 
    << "\n... total # of runs: " << nRunTotal_ 
    << "\n-------  ( of which     " << nNonlTotal_ << " with  NonlinModule.)"
    << "\n... # of commits:    " << nCommTotal_
    << "\n... # of continues:  " << nContTotal_
    << "\n-------  ( of which     " << nCancCont_ << " were canceled.)"
    << "\n... # of cancels:    " << nCancTotal_
    << "\n..."
    << "\n... total # of iterations: " << nIterTotal_
    << "\n-------  ( of which " << nIterCont_ << " after continue.)"
    << endl;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void FlexArclenModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties myProps = props.findProps ( myName_ );

  myProps.find ( doTmpNonlin_, TMP_NONLIN_PROP );
  myProps.find ( writeStats_ , WRITE_STAT_PROP );
  myProps.find ( tsBased_ ,    TIMESTEP_BASED  );
  
  if ( doTmpNonlin_ )
  {
    myProps.find ( keepLargeTmp_, KEEP_LARGE_PROP );
  }

  nonlinSolver_  ->configure ( props, globdat );
  arclenSolver_  ->configure ( props, globdat );
  tsArclenSolver_->configure ( props, globdat );

  if ( writeStats_ )
  {
    cpuTimer_.start ();

    statOut_  = newInstance<PrintWriter>(
                newInstance<FileWriter> ( "flex.stats" ) );
    *statOut_ << "   nRun  nComm  nCanc  nCont  nCaCo  nIter  CPU" << endl;
  }

    
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void FlexArclenModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf  = conf.makeProps ( myName_ );

  myConf.set( "type" , TYPE_NAME );

  myConf.set( TMP_NONLIN_PROP, doTmpNonlin_ );
  myConf.set( WRITE_STAT_PROP, writeStats_  );
  myConf.set( TIMESTEP_BASED,  tsBased_     );


  if ( doTmpNonlin_ )
  {
    myConf.set ( KEEP_LARGE_PROP, keepLargeTmp_ );
  }

  nonlinSolver_  ->getConfig ( conf, globdat );
  arclenSolver_  ->getConfig ( conf, globdat );
  tsArclenSolver_->getConfig ( conf, globdat );
}

//-----------------------------------------------------------------------
//   storeMaxIncrs_
//-----------------------------------------------------------------------

void FlexArclenModule::storeMaxIncrs_

  ( const Properties&       globdat )

{
  // Store maximum values of increment tried in this time step
  // These values are used for increaseStep_

  Properties  params;
  double      inc;

  model_->takeAction ( XActions::GET_STEP_SIZE, params, globdat );

  if ( params.find ( inc, XProps::STEP_SIZE ) )
  {
    if ( currentSolver_ == nonlinSolver_ )
    {
      maxDispTried_ = max ( maxDispTried_, inc );
    }
    else
    {
      maxArcTried_  = max ( maxArcTried_ , inc );
    }
  }
  else
  {
    System::warn() << "Is there a DispArclenModel defined?" << endl;
  }
}

//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------

void FlexArclenModule::commit_

  ( const Properties&       globdat )

{
  // System::out() << "Committing ... doneTmpNonlin " << doneTmpNonlin_
    // << ", nContinues " << nContinues_ << endl;

  // store this solution and proceed with next time step
  // model->takeAction(COMMIT) is called by solver!

  currentSolver_->commit ( globdat );

  ++nCommTotal_;

  // reset variables

  istep0_ = istep_;

  nCancels_       = 0;
  nChanges_       = 0;
  nContinues_     = 0;
  nonlinTriedAll_ = false;
  arclenTriedAll_ = false;
  arclenTriedSma_ = false;
  doneTmpNonlin_  = 0;
  persistent_     = false;
  triedCareful_   = false;
  maxDispTried_   = 0.;
  maxArcTried_    = 0.;

  if ( isTmpNonlin_ )
  {
    // tmpNonlin was succesfull, switch back to Arclen

    switchToArclen_ ( globdat );

    isTmpNonlin_ = false;
  }

  doneResSwitch_ = false;

  if ( isCareful_ )
  {
    Properties params;
    isCareful_ = false;
    model_->takeAction ( XActions::STOP_CAREFUL, params, globdat );
  }

  if ( writeStats_ )
  {
    NumberFormat nf;

    nf.setIntegerWidth   ( 6 );
    nf.setFractionDigits ( 4 );
    nf.setScientific     ( true );

    *statOut_ <<         nf.print ( nRunTotal_  )
              << ' '  << nf.print ( nCommTotal_ )
              << ' '  << nf.print ( nCancTotal_ )
              << ' '  << nf.print ( nContTotal_ )
              << ' '  << nf.print ( nCancCont_  )
              << ' '  << nf.print ( nIterTotal_ ) 
              << "  " << nf.print ( cpuTimer_.toDouble() ) 
              << endl;

    statOut_  -> flush();
  }
}

//-----------------------------------------------------------------------
//   continue_
//-----------------------------------------------------------------------

void FlexArclenModule::continue_

  ( const Properties&       params,
    const Properties&       globdat )

{
  // the solver is not canceled, because the current solution
  // (converged but not accepted) is a good starting point
  
  System::out() << "Continuing with this load step\n\n";

  istep_ = istep0_;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  ++nContinues_;
  ++nContTotal_;

  // find flag that was possibly set in DispArclenModel::checkCommit 
  // (set to true if arcLength_ < swtEnergy_ )

  bool modelAllowsIt = false;

  params.find ( modelAllowsIt, XProps::ALLOW_TMP );

  if ( ( currentSolver_ == arclenSolver_     ||
         currentSolver_ == tsArclenSolver_ ) && 
         doTmpNonlin_                        && 
       ! doneTmpNonlin_                      &&
         modelAllowsIt                    )
  {
    // try switch to displacement control for remainder of time step

    // store current arc-length in variable arcBackup_

    Properties backup;
    model_->takeAction ( XActions::GET_STEP_SIZE, backup, globdat );
    backup.get ( arcBackup_, XProps::STEP_SIZE );

    params.set ( XActions::TEMPORARY, true );

    // allow for switching back based on swtEnergy? 
    // depends on whether you set persistent_ here or not

    // isTmpNonlin_   = persistent_ = switchToNonlin_ ( globdat, params );

    // temporarily reset nonlinTriedAll_

    bool oldVal     = nonlinTriedAll_;
    nonlinTriedAll_ = false;
    isTmpNonlin_    = switchToNonlin_ ( globdat, params );
    nonlinTriedAll_ = oldVal;

    doneTmpNonlin_ = 1;
  }

  // NB: the same params object is used, so ALLOW_TMP is accessible
  // in the model in this action

  params.set ( XProps::N_CONTINUES, nContinues_ );

  model_->takeAction ( XActions::CONTINUE, params, globdat );
}

//-----------------------------------------------------------------------
//   cancel_
//-----------------------------------------------------------------------

void FlexArclenModule::cancel_

  ( const Properties&       globdat )

{
  // discard solution, go back to the beginning of the time step
  // model->takeAction(CANCEL) is called by solver!

  ++nCancels_;
  ++nCancTotal_;

  nCancCont_ += nContinues_;

  // System::out() << " Canceling ... " << doneTmpNonlin_
    // << ", nContinues " << nContinues_ << " -> " << nCancCont_ << endl;

  currentSolver_->cancel ( globdat );

  istep_ = istep0_;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  if ( isCareful_ )
  {
    triedCareful_ = true;
    Properties params;
    isCareful_ = false;
    model_->takeAction ( XActions::STOP_CAREFUL, params, globdat );
  }
  else
  {
    triedCareful_ = false;
  }
}

//-----------------------------------------------------------------------
//   switchSolver_
//-----------------------------------------------------------------------

bool FlexArclenModule::switchSolver_

  ( const Properties&       globdat )

{
  return ( switchToNonlin_(globdat) || switchToArclen_(globdat) );
}

//-----------------------------------------------------------------------
//   switchToNonlin_
//-----------------------------------------------------------------------

bool FlexArclenModule::switchToNonlin_

  ( const Properties&       globdat,
    const Properties&       params )

{
  if ( nonlinTriedAll_ || currentSolver_ == nonlinSolver_ ) return false;

  model_->takeAction ( XActions::TO_DISP, params, globdat );
    
  currentSolver_ = nonlinSolver_;

  // allow for switching back based on swtEnergy? 
  // depends on whether you set persistent_ here or not
      
  if ( arclenTriedSma_ ) persistent_ = true;

  return true;
}

//-----------------------------------------------------------------------
//   switchToArclen_
//-----------------------------------------------------------------------

bool FlexArclenModule::switchToArclen_

  ( const Properties&       globdat,
    const Properties&       params )

{
  if ( arclenTriedAll_ || currentSolver_ == arclenSolver_
       || currentSolver_ == tsArclenSolver_ ) return false;

  model_->takeAction ( XActions::TO_ARCLEN, params, globdat );

  if ( arclenTriedSma_ )
  {
    params.set ( XProps::ADAPT_HOW,  "increase"   );
    params.set ( XProps::ADAPT_FROM, maxArcTried_ );

    model_->takeAction ( XActions::ADAPT_STEP, params, globdat );
  }

  // give loadScale from model to solver 

  //params.get( loadScale, "OldLoadScale" );

  //System::out() << "LoadScale set to " << loadScale << "\n";

  //arclenSolver_->setLoadScale ( loadScale );

  if ( tsBased_ )
  {
    currentSolver_ = tsArclenSolver_;
  }
  else
  {
    currentSolver_ = arclenSolver_;
  }

  return true;
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

bool FlexArclenModule::reduceStep_

  ( const Properties&       globdat )

{
  if ( currentSolver_ == nonlinSolver_ ? 
                         nonlinTriedAll_ :
                         arclenTriedSma_ )
  {
    return false;
  }

  Properties  params;

  bool        reduced = false;

  params.set ( XProps::ADAPT_HOW, "reduce" );

  model_->takeAction ( XActions::ADAPT_STEP, params, globdat );

  params.find ( reduced, XActions::DONE );

  currentSolver_ == nonlinSolver_ ?
                    nonlinTriedAll_ = ! reduced : 
                    arclenTriedSma_ = ! reduced ;

  return reduced;
}

//-----------------------------------------------------------------------
//   increaseStep_
//-----------------------------------------------------------------------

bool FlexArclenModule::increaseStep_

  ( const Properties&       globdat )

{
  Properties  params;
  bool        increased = false;

  if ( ! arclenTriedAll_ )
  {

    // switch solver if not already in Arclen

    switchToArclen_(globdat);

    // increase energy increment

    params.set ( XProps::ADAPT_HOW,  "increase"   );
    params.set ( XProps::ADAPT_FROM, maxArcTried_ );

    model_->takeAction ( XActions::ADAPT_STEP, params, globdat );

    params.find ( increased, XActions::DONE );

    // give it up for the arclength solver if increase fails

    arclenTriedAll_ = ! increased;
  }

  return increased;
}

//-----------------------------------------------------------------------
//   turnOffResSwt_
//-----------------------------------------------------------------------

bool FlexArclenModule::turnOffResSwt_

  ( const Properties&       globdat )

{
  if ( ! doneResSwitch_ ) return false;

  doneResSwitch_  = false;
  arclenTriedAll_ = false;
  noResSwitch_    = true;

  return switchToArclen_ ( globdat );
}

//-----------------------------------------------------------------------
//   considerCarefulMode_
//-----------------------------------------------------------------------

bool FlexArclenModule::considerCarefulMode_

  ( const Properties&       globdat )

{
  Properties actionParams;

  isCareful_ = model_->takeAction 
    ( XActions::BE_CAREFUL, actionParams, globdat );

  triedCareful_ = true;

  return isCareful_;
}

//-----------------------------------------------------------------------
//   getDissForce_
//-----------------------------------------------------------------------

void FlexArclenModule::getDissForce_

  ( const Properties&       globdat ) const

{
  // get fDiss for dissipation based arclength:
  // first initialize, then take action
  // NB: thermal part is added by ThermalModule

  Properties params;
  Vector     fDiss;

  globdat.find ( fDiss, XProps::DISSIPATION_FORCE );

  Ref<DofSpace> dofs = DofSpace::find ( globdat );

  if ( fDiss.size() != dofs->dofCount() )
  {
    fDiss.resize ( dofs->dofCount() );
    globdat.set  ( XProps::DISSIPATION_FORCE, fDiss );
  }

  fDiss = 0.;

  model_->takeAction ( XActions::GET_DISS_FORCE, params, globdat );
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  FlexArclenModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat )

{
  Properties  myConf  = conf .makeProps ( name );
  Properties  myProps = props.findProps ( name );

  Ref<NonlinModule>   solver1 = newInstance<NonlinModule> 

                                  ( joinNames ( name, NONLIN ) );

  Ref<ArclenModule>   solver2 = newInstance<ArclenModule> 

                                  ( joinNames ( name, ARCLEN ) );

  Ref<TSArclenModule> solver3 = newInstance<TSArclenModule> 

                                  ( joinNames ( name, TSARCLEN ) );

  return newInstance<Self> ( name, solver1, solver2, solver3 );
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareFlexArclenModule
//-----------------------------------------------------------------------

void declareFlexArclenModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( FlexArclenModule::TYPE_NAME,
                         & FlexArclenModule::makeNew );
}

