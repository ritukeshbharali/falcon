
/** @file ReduceSteppingModule.cpp
 *  @brief Implements (time) stepping module with decreasing step-size.
 *  
 *  This class implements a module for reducing (time) 
 *  step sizes based on user input or when the solver
 *  fails to converge. Unlike the Adaptive Stepping 
 *  Module, the step size is never increased based on
 *  optimal iterations per step. Once the stepsize is 
 *  changed, the module throws the action 
 *  XActions::SET_STEP_SIZE. Models, such as the 
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


/* Include jem and jive headers */

#include <jem/base/limits.h>
#include <jem/base/System.h>
#include <jem/base/Exception.h>
#include <jem/base/IllegalOperationException.h>
#include <jem/io/FileWriter.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/util/Globdat.h>
#include <jive/util/utilities.h>
#include <jive/util/DofSpace.h>
#include <jive/implict/SolverInfo.h>
#include <jive/implict/SolverModule.h>
#include <jive/app/ModuleFactory.h>
#include <jive/model/Actions.h>

/* Include user-defined class headers */

#include "ReduceSteppingModule.h"
#include "util/XNames.h"


using jem::IllegalOperationException;
using jem::System;
using jem::max;
using jem::maxOf;
using jem::newInstance;
using jem::io::endl;
using jem::io::FileWriter;
using jem::util::StringUtils;
using jive::Vector;
using jive::implict::SolverInfo;
using jive::util::Globdat;
using jive::util::joinNames;
using jive::util::DofSpace;
using jive::implict::newSolverModule;


//=======================================================================
//   class ReduceSteppingModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* ReduceSteppingModule::TYPE_NAME        = "ReduceStepping";
const char* ReduceSteppingModule::SOLVER           = "solver";
const char* ReduceSteppingModule::START_INCR_PROP  = "startIncr";
const char* ReduceSteppingModule::MIN_INCR_PROP    = "minIncr";
const char* ReduceSteppingModule::REDUCE_STEP_PROP = "reduceStep";
const char* ReduceSteppingModule::REDUCTION_PROP   = "reduction";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

ReduceSteppingModule::ReduceSteppingModule 

  ( const String&  name,
    Ref<SolverModule>  solver ) :

      Super       ( name   ),
      solver_     ( solver ),
      istep_      ( 0 ),
      istep0_     ( 0 ),
      startIncr_  ( 1. ),
      minIncr_    ( 1.e-5 ),
      cutStep_    ( -1 ),
      reduction_  ( 0.45 ),
      increment0_ ( 0. ),
      increment_  ( 0. )

{}

ReduceSteppingModule::~ReduceSteppingModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status ReduceSteppingModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  model_ = Model::get ( globdat, getContext() );
  istep_ = istep0_ = 0;

  globdat.set ( Globdat::TIME_STEP, istep_ );
  globdat.set ( "var.accepted", false      );

  solver_->init ( conf, props, globdat );

  return OK;
}


//-----------------------------------------------------------------------
//   run
//-----------------------------------------------------------------------


Module::Status ReduceSteppingModule::run ( const Properties& globdat )

{
  using jem::Exception;

  Properties  info      = SolverInfo::get ( globdat );
  Properties  params;

  bool        accept    = true;
  bool        discard   = false;
  bool        converged = false;
  String      terminate;

  info.clear ();

  setStepSize_( globdat );

  solver_->advance ( globdat );
  
  // Advance the time/load step number

  istep_ = istep0_ + 1;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  // try to find a solution

  try
  {
    solver_->solve ( info, globdat );

    converged = true;
  }
  catch ( const jem::Exception &ex )
  {
    System::out() << "ReduceSteppingModule caught a solver exception: \n" 
                  << "   " << ex.where() << "\n"
                  << "   " << ex.what() << "\n\n";
                  
    if ( ex.what() == "invalid matrix element(s): NaN" )
    {
      //return EXIT;
    }

    discard    = true;
  }

  if ( converged )
  {
    // Ask the model whether to accept converged solution
    // (Model can ask for cancel via DISCARD)

    model_->takeAction ( XActions::CHECK_COMMIT, params, globdat );

    params.find ( discard, XActions::DISCARD );
    params.find ( accept,  XActions::ACCEPT  );
  }

  // set flag in globdat for output indicating solution quality

  globdat.set ( "var.accepted", ( accept && ! discard ) );

  if ( ! discard )
  {
    // take measures based on result of CHECK_COMMIT

    if ( accept )
    {
      // The final solution for this time step has been obtained

      commit_ ( globdat );

      if ( params.find ( terminate, XActions::TERMINATE ) ) 
      {
        if ( terminate == "sure" )
        {
          System::out() << "Some model gave the terminate signal during "
            << "checkCommit,\nReduceSteppingModule terminates run.\n"
            << "NB: output modules are not called for this time step\n";
          return EXIT;
        }
      }
    }
    else
    {
      // We have an equilibrium solution but it is not accepted:
      // continue with this time step

      continue_ ( params, globdat );
    }
  }
  else
  {
    // 1. discard solution

    cancel_ ( globdat );

    // 2. find proper strategy to try again
   
    // try to reduce the step size

    if ( reduceStep_ ( globdat ) )
    {
      System::out() << "Returning to the same load step "
                    << "with reduced step\n";
    }

    // give it up

    else
    {
      System::out() << "\n*** ReduceSteppingModule is out of inspiration,\n"
                    << "    can't find a winning strategy for this step."
                    << " Sorry.\n\n";

      return EXIT;
    }
  }
  
  
  // check if the minimum increment is reached

  if ( increment_ < minIncr_ )
  {
    System::out() << "\n=======================================================\n";
    System::out() << "ReduceSteppingModule reached minimum increment.\nSimulation will now be terminated.\n";
    System::out() << "=======================================================\n\n";
    return EXIT;
  }
  
  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void ReduceSteppingModule::shutdown ( const Properties& globdat )
{
  solver_->shutdown ( globdat );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ReduceSteppingModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties myProps = props.findProps ( myName_ );
  
  solver_->configure ( props, globdat );

  myProps.find ( startIncr_,   START_INCR_PROP  );
  myProps.find ( minIncr_,     MIN_INCR_PROP    );
  myProps.find ( cutStep_ , REDUCE_STEP_PROP );
  myProps.find ( reduction_,   REDUCTION_PROP   );

  increment_ = startIncr_;

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ReduceSteppingModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf  = conf.makeProps ( myName_ );

  myConf.set ( "type", TYPE_NAME );
  
  myConf.set ( START_INCR_PROP,  startIncr_  );
  myConf.set ( MIN_INCR_PROP  ,  minIncr_    );
  myConf.set ( REDUCE_STEP_PROP, cutStep_ );
  myConf.set ( REDUCTION_PROP ,  reduction_  );

  solver_->getConfig ( conf, globdat );
}

//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------

void ReduceSteppingModule::commit_

  ( const Properties&       globdat )

{
  // store this solution and proceed with next time step
  // model->takeAction(COMMIT) is called by solver!

  solver_->commit ( globdat );

  // update old step and old step increment

  istep0_     = istep_;
  increment0_ = increment_;

  // reduce step if a certain step is reached

  if ( istep_ == cutStep_ )
  {
    increment_ *= reduction_;
  }

}

//-----------------------------------------------------------------------
//   continue_
//-----------------------------------------------------------------------

void ReduceSteppingModule::continue_

  ( const Properties&       params,
    const Properties&       globdat )

{
  // the solver is not canceled, because the current solution
  // (converged but not accepted) is a good starting point
  
  System::out() << "Continuing with this load step\n\n";

  istep_ = istep0_;

  globdat.set ( Globdat::TIME_STEP, istep_ );

  model_->takeAction ( XActions::CONTINUE, params, globdat );
}

//-----------------------------------------------------------------------
//   cancel_
//-----------------------------------------------------------------------

void ReduceSteppingModule::cancel_

  ( const Properties&       globdat )

{
  // discard solution, go back to the beginning of the time step
  // model->takeAction(CANCEL) is called by solver!

  solver_->cancel ( globdat );

  istep_ = istep0_;

  globdat.set ( Globdat::TIME_STEP, istep_ );
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

bool ReduceSteppingModule::reduceStep_

  ( const Properties&       globdat )

{
  increment_ *= reduction_;

  return true;
}

//-----------------------------------------------------------------------
//   setStepSize()
//-----------------------------------------------------------------------

void ReduceSteppingModule::setStepSize_

  ( const Properties& globdat ) const

{
  using jive::model::ActionParams;
  using jive::model::Actions;
  
  Properties params;

  params.set ( XProps::STEP_SIZE  , increment_  );
  params.set ( XProps::STEP_SIZE_0, increment0_ );

  // set step size in models
  model_->takeAction ( XActions::SET_STEP_SIZE, params, globdat );
  
  // set step size in solverModule
  Properties props;
  props.set ( joinNames ( myName_ , ".solver.deltaTime" ) , increment_ ); 
  solver_->configure ( props, globdat );
  
}

//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------

Ref<Module>  ReduceSteppingModule::makeNew

  ( const String&           name,
    const Properties&       conf,
    const Properties&       props,
    const Properties&       globdat )

{
  Properties  myConf  = conf .makeProps ( name );
  Properties  myProps = props.findProps ( name );

  Ref<SolverModule> solver = newSolverModule 
    ( joinNames ( name, SOLVER ), conf, props, globdat  );

  return newInstance<Self> ( name, solver );
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareReduceSteppingModule
//-----------------------------------------------------------------------

void declareReduceSteppingModule ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( ReduceSteppingModule::TYPE_NAME,
                         & ReduceSteppingModule::makeNew );
}