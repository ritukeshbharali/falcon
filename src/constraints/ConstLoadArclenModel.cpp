/*
 * 
 *  Copyright (C) 2009 TU Delft. All rights reserved.
 *  
 *  This class implements a model that combines displacement and 
 *  arclen control. This model is designed for use in combination 
 *  with FlexArclenModule to have a flexible path following solver. 
 *  Boundary conditions are applied via a DispArclenBCs object that 
 *  is member of this class.
 *
 *  There are two different solution strategies, each with it's own 
 *  solver (NonlinModule or ArclenModule) and way of handling 
 *  boundary conditions. The FlexArclenModule decides which strategy 
 *  is applied and communicates with this model via takeAction calls. 
 *  
 *  The tasks of this model are:
 *   - let the DispArclenBCs apply correct boundary conditions
 *     (adaptive in step type and step size)
 *   - perform actions required for ArclenModule related to energy
 *     arclength method (compute dissipation etc)
 *   - pass information to FlexArclenModule for strategy decision 
 *     making
 *
 *  Authors: 
 *  V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  F.P. van der Meer, F.P.vanderMeer@tudelft.nl
 *         
 *  Date: 27 February 2009
 *
 */

// TO-DO: do not use the unit load from this model if
//        other (FE) model has updated it. implement
//        just like bool assem. might need more XNames.


#include <jem/base/array/operators.h>
#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/utilities.h>
#include <jem/util/Event.h>
#include <jem/util/StringUtils.h>
#include <jem/util/ArrayBuffer.h>
#include <jive/util/error.h>
#include <jive/util/Printer.h>
#include <jive/util/utilities.h>
#include <jive/util/ItemSet.h>
#include <jive/util/DofSpace.h>
#include <jive/util/Constraints.h>
#include <jive/util/VectorManager.h>
#include <jive/util/Globdat.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/implict/ArclenActions.h>
#include <jive/implict/SolverInfo.h>

#include "ConstLoadArclenModel.h"
#include "util/XNames.h"


extern "C"
{
  #include <math.h>
}

using jem::io::endl;
using jem::numeric::axpy;
using jem::util::StringUtils;
using jem::util::ArrayBuffer;
using jive::IdxVector;
using jive::util::VectorManager;
using jive::model::Actions;
using jive::model::ActionParams;
using jive::model::StateVector;
using jive::implict::ArclenActions;
using jive::implict::ArclenParams;
using jive::implict::SolverInfo;


//=======================================================================
//   class ConstLoadArclenModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  ConstLoadArclenModel::TYPE_NAME       = "ConstLoadArclen";

const char*  ConstLoadArclenModel::OPT_ITER_PROP   = "optIter";
const char*  ConstLoadArclenModel::SWT_ITER_PROP   = "swtIter";
const char*  ConstLoadArclenModel::SWT_ENER_PROP   = "swtEnergy";
const char*  ConstLoadArclenModel::MIN_INCR_PROP   = "minIncr";
const char*  ConstLoadArclenModel::MAX_INCR_PROP   = "maxIncr";
const char*  ConstLoadArclenModel::REDUCTION_PROP  = "reduction";

const char*  ConstLoadArclenModel::DELTA_STATE_    = "deltaState0";

const idx_t  ConstLoadArclenModel::U_LOAD_         = 1 << 0;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


ConstLoadArclenModel::ConstLoadArclenModel

  ( const String&      name,
    const Ref<Model>&  child ) :

    Super   ( name  ),
    out_    ( System::info( name ) )

{
  updated_       = 0;
  maxNIter_      = 0;
  optIter_       = 4;
  minIncr_       = 1.0e-3;
  maxIncr_       = 1.0e+1;
  arcLength_     = 0.0;
  arcLength0_    = 0.0;
  lastArclength_ = 0.0;
  reduction_     = .55;
  swtEner_       = Float::MAX_VALUE;
  swtIter_       = 100;

  dtime_         = 1.0;
  dtime0_        = 1.0;

  gamma_         = 1.e+20;

  isLoadControl_ = true; 
  isTmpLoad_     = false; 
  onceDown_      = false;
  hasDespaired_  = false;
}


ConstLoadArclenModel::~ConstLoadArclenModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool ConstLoadArclenModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  // initialization

  if ( action == Actions::INIT )
  {
    init_ ( globdat );

    return true;
  }

  // compute arclength function

  if ( action == ArclenActions::GET_ARC_FUNC )
  {
    evalArcFunc_ ( params, globdat );

    return true;
  }

  if ( action == ArclenActions::GET_UNIT_LOAD )
  {
    getUnitLoad_ ( params, globdat );

    return true;
  }

  // perform checks on converged solution

  if ( action == XActions::CHECK_COMMIT )
  {
    checkCommit_ ( params, globdat );

    checkSwitch_ ( params, globdat );

    reportProgress_();

    return true;
  }

  // proceed to next time step

  if ( action == Actions::COMMIT )
  {
    commit_ ( params, globdat );

    return true;
  }

  // reset variables

  if ( action == Actions::CANCEL )
  {
    cancel_ ( params, globdat );

    return true;
  }

  // set step size (dEnergy or du)

  if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params, globdat );

    return true;
  }

  // adapt step size

  if ( action == XActions::ADAPT_STEP )
  {
    String how;

    params.get ( how, XProps::ADAPT_HOW );

    if      ( how == "reduce" )
    {
      reduceStep_ ( params, globdat );
    }
    else if ( how == "increase" )
    {
      increaseStep_ ( params, globdat  );
    }

    return true;
  }

  // switch to arclength control

  if ( action == XActions::TO_ARCLEN )
  {
    toArclControl_ ( params, globdat );

    return true;
  }

  // switch to displacement control
  
  if ( action == XActions::TO_DISP )
  {
    toLoadControl_ ( params, globdat );

    return true;
  }

  // STOP_ARCL is called by XArclenModule, return DONE to make known 
  // that this is a model that allows for strategy switching
  // Switch will be issued by FlexArclenModule

  if ( action == XActions::STOP_ARCLEN )
  {
    params.set ( XActions::DONE, true );

    return true;
  }

  if ( action == "DESPAIR" ) 
  {
    hasDespaired_ = true;

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   findModel
//-----------------------------------------------------------------------


Model* ConstLoadArclenModel::findModel ( const String& name ) const
{
  return const_cast<Self*> ( this );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void ConstLoadArclenModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  double maxD = Float::MAX_VALUE;

  myProps.find ( optIter_,   OPT_ITER_PROP,  1,        1000 );
  myProps.find ( swtIter_,   SWT_ITER_PROP,  1,        1000 );
  myProps.find ( swtEner_,   SWT_ENER_PROP,  0.0,      maxD );
  myProps.find ( minIncr_,   MIN_INCR_PROP,  0.0,      maxD );
  myProps.find ( maxIncr_,   MAX_INCR_PROP,  minIncr_, maxD );
  myProps.find ( reduction_, REDUCTION_PROP, 0.0,      1.0  );

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void ConstLoadArclenModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( OPT_ITER_PROP,    optIter_      );
  myConf.set ( SWT_ITER_PROP,    swtIter_      );
  myConf.set ( SWT_ENER_PROP,    swtEner_      );
  myConf.set ( MIN_INCR_PROP,    minIncr_      );
  myConf.set ( MAX_INCR_PROP,    maxIncr_      );
  myConf.set ( REDUCTION_PROP,   reduction_    );
}


//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void ConstLoadArclenModel::setMaxIter ( idx_t count )
{
  JEM_PRECHECK ( count > 0 );

  optIter_ = count;
}


//-----------------------------------------------------------------------
//   setIncrRange
//-----------------------------------------------------------------------


void ConstLoadArclenModel::setIncrRange

  ( double  minIncr,
    double  maxIncr )

{
  JEM_PRECHECK ( minIncr >= 0.0 && minIncr < maxIncr );

  minIncr_ = minIncr;
  maxIncr_ = maxIncr;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Model> ConstLoadArclenModel::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return newInstance<Self> ( name );
}


//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::init_ ( const Properties& globdat )
{
  using jem::util::connect;

  dofs_    = DofSpace    :: get ( globdat, getContext() );
  updated_ = 0;

  connect ( dofs_->newSizeEvent,   this, & Self::dofsChanged_ );
  connect ( dofs_->newOrderEvent,  this, & Self::dofsChanged_ );

  dofs_->resetEvents ();
}


//-----------------------------------------------------------------------
//   initLoad_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::initLoad_ ( const Properties& globdat )
{
  // In the current implementation, it does nothing. This is because
  // the unitLoad is computed by the FE model.
}


//-----------------------------------------------------------------------
//   evalArcFunc_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::evalArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // In this implementation, arc-length function is computed by getting
  // the dissipation computed from the FE model and subtracting the
  // arc-length value.

  double  fvalue;
  double  dEnergy;
  double  jac11;

  bool fromFEModel = globdat.find ( dEnergy, XProps::FE_DISSIPATION );

  if ( !fromFEModel )
  {
    throw IllegalInputException ( JEM_FUNC,
          "FE Model must update FE_DISSIPATION!" );
  }

  // ( Original implementation )

  fvalue = dEnergy - arcLength_;

  params.set ( ArclenParams::ARC_FUNC,   fvalue );
}


//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::getUnitLoad_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // In this implementation, unitLoad is taken from the FE Model and
  // set.

  using jive::util::sizeError;

  Vector f, unitLoad, fext;

  bool fromFEModel = globdat.find ( f, XProps::FE_UNIT_LOAD );

  if ( !fromFEModel )
  {
    throw IllegalInputException ( JEM_FUNC,
          "FE Model must update UNIT_LOAD!" );
  }

  params.get ( unitLoad, ArclenParams::UNIT_LOAD );

  if ( f.size() != unitLoad.size() )
  {
    sizeError ( getContext(),
                "unit load vector", f.size(), unitLoad.size() );
  }

  // globdat.get( fext, "fexternal" );
  
  // ( When pressure equation is multiplied by dtime )
  // unitLoad = -f + fext;

  // ( When pressure equation is not multiplied by dtime )
  unitLoad = -f; // Works for hydraulic fracture problems
  //unitLoad = fext;
}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{

  // If load control method is adopted, do nothing. This means that we
  // continue with the same step-size and the constant external load.
  // ( TO-DO: Add provision to change step-size! )

  // If Arclen control method is adopted, depending on the iterations,
  // change the arc-length value.

  if ( !isLoadControl_ )
  {

    arcLength0_ = arcLength_;
    dtime0_     = dtime_;

    if ( maxNIter_ < optIter_ )
    {
      arcLength_ *=  1.15; // 1.25;
    }
    else
    {
      double n    = ( maxNIter_ - optIter_ ) / 4.0;
      arcLength_ *= ::pow ( 0.5, n );
    }

    // bound the arc-length value

    if        ( arcLength_ < minIncr_ )
    {
      arcLength_ = minIncr_;
    }
    else if   ( arcLength_ > maxIncr_ )
    {
      arcLength_ = maxIncr_;
    }

    out_ << "arcLength = " << arcLength_ << endl;

  }

  maxNIter_     = 0;
  hasDespaired_ = false;
  tempStep_     = false;
}

//-----------------------------------------------------------------------
//   cancel_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::cancel_

  ( const Properties&  params,
    const Properties&  globdat )

{
  maxNIter_     = 0;
  hasDespaired_ = false;
  dtime_        = dtime0_;

  if ( arcLength_ < 1.e-2 * minIncr_ )
  {
    params.set ( XActions::TERMINATE, "sure" );
  }
}

//-----------------------------------------------------------------------
//   reportProgress_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::reportProgress_ () const

{
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::checkCommit_

  ( const Properties&  params,
    const Properties&  globdat )

{
}

//-----------------------------------------------------------------------
//   checkSwitch_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::checkSwitch_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Globdat;

  // store number of iterations for adaptive stepping

  idx_t  iiter;

  params.get ( iiter, SolverInfo::ITER_COUNT );
  
  int           step;
  globdat.get ( step, Globdat::TIME_STEP );

  maxNIter_ = max ( maxNIter_, iiter );

  if ( isLoadControl_ && ! tempStep_ && step > 1 )
  {
    // get the dissipated energy in this step
    // FE_Dissipation must be updated whenever GET_MATRIX0 is called
    // in the FE model.

    bool fromFEModel = globdat.find ( lastArclength_, XProps::FE_DISSIPATION );

    if ( !fromFEModel )
    {
      throw IllegalInputException ( JEM_FUNC,
            "FE Model must update FE_DISSIPATION!" );
    }

    if ( lastArclength_ > swtEner_ ||
        ( iiter > swtIter_ && lastArclength_ > minIncr_ ) ) 
    {
      params.set ( XActions::DO_SWITCH, "please" );

      arcLength_ = lastArclength_;
    }
  }
  else
  {
    // TO-DO: Explore options to switch back to load control!
  }


}

//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::setStepSize_

  ( const Properties&  params,
    const Properties&  globdat )

{
  double dt;
  params.get ( dt, XProps::STEP_SIZE   );

  dtime_ = dt;
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::reduceStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( !isLoadControl_ )
  {
    arcLength_ *= reduction_;

    out_ << "Path following parameter reduced to "
         << arcLength_ << "\n";

    params.set ( XActions::DONE, true );
  }
}

//-----------------------------------------------------------------------
//   increaseStep_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::increaseStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( !isLoadControl_ )
  {
    arcLength_ /= reduction_;

    out_ << "Path following parameter increased to "
         << arcLength_ << "\n";

    params.set ( XActions::DONE, true );
  }
}

//-----------------------------------------------------------------------
//   toArclControl_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::toArclControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  isTmpLoad_ = isLoadControl_ = false;

  out_ << "Switching to arclength control with dEnergy: " 
       << arcLength_ << endl;
}

//-----------------------------------------------------------------------
//   toLoadControl_
//-----------------------------------------------------------------------

void ConstLoadArclenModel::toLoadControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  isLoadControl_ = true;
  out_ << "Switching back to load control " << "\n";
}


//-----------------------------------------------------------------------
//   dofsChanged_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::dofsChanged_ ()
{
  updated_ = 0;
}


//-----------------------------------------------------------------------
//   consChanged_
//-----------------------------------------------------------------------


void ConstLoadArclenModel::consChanged_ ()
{
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareConstLoadArclenModel
//-----------------------------------------------------------------------


void declareConstLoadArclenModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( ConstLoadArclenModel::TYPE_NAME,
                          & ConstLoadArclenModel::makeNew );
}
