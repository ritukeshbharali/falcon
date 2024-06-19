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

#include "DispArclenModel.h"
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
//   class DispArclenModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  DispArclenModel::TYPE_NAME       = "DispArclen";

const char*  DispArclenModel::CONSTRAINT_PROP = "constraints";
const char*  DispArclenModel::OPT_ITER_PROP   = "optIter";
const char*  DispArclenModel::SWT_ITER_PROP   = "swtIter";
const char*  DispArclenModel::SWT_ENER_PROP   = "swtEnergy";
const char*  DispArclenModel::MIN_INCR_PROP   = "minIncr";
const char*  DispArclenModel::MAX_INCR_PROP   = "maxIncr";
const char*  DispArclenModel::DISP_INCR_PROP  = "dispIncr";
const char*  DispArclenModel::INIT_DISP_PROP  = "initDisp";
const char*  DispArclenModel::MIN_DISP_PROP   = "minDispIncr";
const char*  DispArclenModel::MAX_DISP_PROP   = "maxDisp";
const char*  DispArclenModel::MAX_DISS_PROP   = "maxTotalDiss";
const char*  DispArclenModel::REDUCTION_PROP  = "reduction";
const char*  DispArclenModel::PREFER_DC_PROP  = "preferDisp";

const char*  DispArclenModel::DELTA_STATE_    = "deltaState0";

const idx_t  DispArclenModel::U_LOAD_         = 1 << 0;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


DispArclenModel::DispArclenModel

  ( const String&      name,
    const Ref<Model>&  child ) :

    Super   ( name  ),
    dispBC_ ( newInstance<DispArclenBCs>() ),
    out_    ( System::info( name ) )

{
  updated_       = 0;
  maxNIter_      = 0;
  optIter_       = 4;
  minIncr_       = 1.0e-3;
  maxIncr_       = 1.0e+1;
  dispIncr0_     = 1.0;
  arcLength_     = 0.0;
  lastArcl_      = 0.0;
  reduction_     = .55;
  swtEner_       = Float::MAX_VALUE;
  swtIter_       = 100;
  initDisp_      = 0.0;

  totalDiss_     = 0.0;
  minDispIncr_   = 0.0;
  maxDispVal_    = Float::MAX_VALUE;
  maxTotalDiss_  = Float::MAX_VALUE;

  isDispControl_ = true; 
  isTmpDisp_     = false; 
  onceDown_      = false;
  preferDisp_    = false;
  hasDespaired_  = false;
}


DispArclenModel::~DispArclenModel ()
{}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DispArclenModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  // initialization

  if ( action == Actions::INIT )
  {
    dispBC_->initConstraints( globdat );

    dispBC_->setInitVal ( initDisp_ );

    init_ ( globdat );

    return true;
  }

  // store the load scale whenever it has been updated by other models
  // [during displacement control, for possible future switch]

  if ( action == Actions::GET_MATRIX0 ||
       action == Actions::GET_INT_VECTOR )
  {
    Vector       fint;
    params.get ( fint, ActionParams::INT_VECTOR );

    dispBC_->storeLoadScale ( globdat, fint );

    return true;
  }

  // compute the external force vector
  // (only necessary with load control)

  if ( action == Actions::GET_EXT_VECTOR  )
  {
    getExtVector_ ( params, globdat );

    return true;
  }

  // apply displacement increment

  if ( action == Actions::GET_CONSTRAINTS )
  {
    if ( isDispControl_ )
    {
      dispBC_->applyConstraints ( params, globdat );
    }

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

    dispBC_->commit( globdat );

    return true;
  }

  // advance to next time step

  if ( action == Actions::ADVANCE )
  {
    if ( isDispControl_ )
    {
      dispBC_->advance( globdat );
    }

    return true;
  }

  // reset variables

  if ( action == Actions::CANCEL )
  {
    cancel_ ( params, globdat );

    return true;
  }

  // get step size (dEnergy or du)

  if ( action == XActions::GET_STEP_SIZE )
  {
    getStepSize_ ( params, globdat );

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

    dispBC_->toArclControl ( params, globdat );

    return true;
  }

  // switch to displacement control
  
  if ( action == XActions::TO_DISP )
  {
    toDispControl_ ( params, globdat );

    dispBC_->toDispControl ( params, globdat );

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

  // GET_DISSIPATION is used for output

  if ( action == "GET_DISSIPATION" )
  {
    StringVector names ( StringUtils::split( myName_, '.' ) );
    params.set ( names.back(), totalDiss_ );

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


Model* DispArclenModel::findModel ( const String& name ) const
{
  return const_cast<Self*> ( this );
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DispArclenModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps = props.findProps ( myName_ );

  double maxD = Float::MAX_VALUE;

  myProps.find ( optIter_, OPT_ITER_PROP, 1,        1000 );
  myProps.find ( swtIter_, SWT_ITER_PROP, 1,        1000 );
  myProps.find ( swtEner_, SWT_ENER_PROP, 0.0,      maxD );
  myProps.find ( minIncr_, MIN_INCR_PROP, 0.0,      maxD );
  myProps.find ( maxIncr_, MAX_INCR_PROP, minIncr_, maxD );

  myProps.find ( reduction_, REDUCTION_PROP,   0.0, 1.0 );
  myProps.find ( preferDisp_, PREFER_DC_PROP );

  Properties childProps = myProps.findProps ( CONSTRAINT_PROP );

  dispBC_->configure ( childProps, globdat );

  myProps.get  ( dispIncr0_, DISP_INCR_PROP );
  myProps.find ( initDisp_,  INIT_DISP_PROP );

  dispBC_->setDispIncr ( dispIncr0_ );

  minDispIncr_ = jem::numeric::abs ( dispIncr0_ );

  myProps.find ( minDispIncr_,   MIN_DISP_PROP, 0.0, minDispIncr_    );
  myProps.find ( maxDispVal_,    MAX_DISP_PROP, 0.0, maxD            );
  myProps.find ( maxTotalDiss_,  MAX_DISS_PROP, 0.0, maxD            );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DispArclenModel::getConfig

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

  myConf.set ( DISP_INCR_PROP, dispIncr0_    );
  myConf.set ( INIT_DISP_PROP, initDisp_     );
  myConf.set ( MIN_DISP_PROP,  minDispIncr_  );
  myConf.set ( MAX_DISP_PROP,  maxDispVal_   );
  myConf.set ( MAX_DISS_PROP,  maxTotalDiss_ );

  myConf.set ( PREFER_DC_PROP,   preferDisp_   );

  Properties childConf = myConf.makeProps ( CONSTRAINT_PROP );
  dispBC_->getConfig ( childConf, globdat );
}


//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void DispArclenModel::setMaxIter ( idx_t count )
{
  JEM_PRECHECK ( count > 0 );

  optIter_ = count;
}


//-----------------------------------------------------------------------
//   setIncrRange
//-----------------------------------------------------------------------


void DispArclenModel::setIncrRange

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


Ref<Model> DispArclenModel::makeNew

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


void DispArclenModel::init_ ( const Properties& globdat )
{
  dofs_    = DofSpace    :: get ( globdat, getContext() );
  updated_ = 0;

  connect_ ();
}


//-----------------------------------------------------------------------
//   initLoad_
//-----------------------------------------------------------------------


void DispArclenModel::initLoad_ ( const Properties& globdat )
{
  // In the current implementation, this is absolutely ridiculous:
  // a vector with size dofCount and all 0.0's except one entry with 1.0
  // and then in getReleasedEnergy, this one entry is selected 
  // and used for multiplication
  // 
  // However, if you want to generalize to more complex load vectors,
  // this is the way to go

  Properties        params;
  ArrayBuffer<idx_t>  idofbuf;

  idx_t ndofs = dofs_->dofCount();
  
  unitLoad_.resize ( ndofs );

  unitLoad_ = 0.0;

  params.set ( ActionParams::EXT_VECTOR, unitLoad_ );

  dispBC_->getUnitLoad ( params, globdat );

  for ( idx_t i = 0; i < ndofs; ++i )
  {
    if ( jem::numeric::abs ( unitLoad_[i] ) > 0. ) 
    {
      idofbuf.pushBack ( i );
    }
  }

  idofs_.ref ( idofbuf.toArray() );

  dofs_ ->resetEvents ();

  updated_ |= U_LOAD_;
}


//-----------------------------------------------------------------------
//   evalArcFunc_
//-----------------------------------------------------------------------


void DispArclenModel::evalArcFunc_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jem::max;
  using jem::isTiny;

  double  fvalue;
  double  dEnergy;

  if ( vspace_ == NIL )
  {
    vspace_ = VectorSpace::get ( dofs_, globdat );
  }

  bool assem = globdat.find ( dEnergy, XProps::FE_DISSIPATION );

  if ( ! assem )
  {    
    dEnergy = getReleasedEnergy_ ( params, globdat );
  }

  fvalue = dEnergy - arcLength_ ;

  params.set ( ArclenParams::ARC_FUNC,   fvalue );

  if ( abs(dEnergy) < 1.e-10 )
  {
    params.set ( ArclenParams::JACOBIAN11, 1.e-4 );
  }
}


//-----------------------------------------------------------------------
//   getUnitLoad_
//-----------------------------------------------------------------------


void DispArclenModel::getUnitLoad_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::sizeError;

  Vector  f;

  bool assem = globdat.find ( f, XProps::FE_UNIT_LOAD );

  if ( ! assem )
  {

    if ( ! (updated_ & U_LOAD_) )
    {
      initLoad_ ( globdat );
    }

    params.get ( f, ArclenParams::UNIT_LOAD );

    if ( f.size() != unitLoad_.size() )
    {
      sizeError ( getContext(),
                  "unit load vector", f.size(), unitLoad_.size() );
    }

    f = unitLoad_;
  }
}

//-----------------------------------------------------------------------
//   getExtVector_
//-----------------------------------------------------------------------


void DispArclenModel::getExtVector_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::sizeError;

  Vector  f;

  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_ ( globdat );
  }

  params.get ( f, ActionParams::EXT_VECTOR );

  if ( f.size() != unitLoad_.size() )
  {
    sizeError ( getContext(),
                "unit load vector", f.size(), unitLoad_.size() );
  }

  f = dispBC_->getLoadScale() * unitLoad_;
}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void DispArclenModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( isDispControl_ )
  {
    // store actual dissipation for possible switch 

    arcLength_ = lastArcl_;
  }

  // store total dissipated energy (for postprocessing)

  totalDiss_ += lastArcl_;

  if ( isDispControl_ && ! isTmpDisp_ && ! hasDespaired_ )
  {
    // adjust increment according to number of iterations

    double oldIncr = dispBC_->getDispIncr();
    double n       = ( maxNIter_ - optIter_ ) / 4.0;
    double newIncr = oldIncr * ::pow ( 0.5, n );

    // check if it falls within allowable interval (defined by user)

    if ( newIncr / dispIncr0_ > 1. )
    {
      newIncr = dispIncr0_;
    }
    else if ( jem::numeric::abs(newIncr) < minDispIncr_ )
    {
      newIncr = dispIncr0_;
    }
    dispBC_->setDispIncr ( newIncr );
  }
  else if ( ! hasDespaired_ )
  {

    if ( maxNIter_ < optIter_ )
    {
      arcLength_ *= 1.2;
    }
    else
    {
      // adjust increment according to number of iterations

      double n    = ( maxNIter_ - optIter_ ) / 4.0;
      arcLength_ *= ::pow ( 0.5, n );

    }

    // check if it falls within allowable interval (defined by user)

    if        ( arcLength_ < minIncr_ )
    {
      arcLength_ = minIncr_;
    }
    else if   ( arcLength_ > maxIncr_ )
    {
      arcLength_ = maxIncr_;
    }
  }
  maxNIter_ = 0;

  hasDespaired_ = false;

  tempStep_ = false;
}

//-----------------------------------------------------------------------
//   cancel_
//-----------------------------------------------------------------------

void DispArclenModel::cancel_

  ( const Properties&  params,
    const Properties&  globdat )

{
  maxNIter_ = 0;
  hasDespaired_ = false;
}

//-----------------------------------------------------------------------
//   reportProgress_
//-----------------------------------------------------------------------

void DispArclenModel::reportProgress_ () const

{
  if ( isDispControl_ )
  {
    out_ << "DispControl: ";
  }
  else
  {
    out_ << "EnergyArclen: ";
  }
  out_ << "dEnergy "     << lastArcl_               << 
          ", loadScale " << dispBC_->getLoadScale() << 
          ", dispVal "   << dispBC_->getDispValue() << endl;
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void DispArclenModel::checkCommit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // terminate the computation if displacement exceeds maximum.
  // be careful with this!

  if ( jem::numeric::abs ( dispBC_->getDispValue() ) > maxDispVal_ ) 
  {
    System::out() << myName_ << " says: TERMINATE because "
      << " disp > maxDispVal." << endl;
    params.set ( XActions::TERMINATE, "sure" );
  }
  else if ( totalDiss_ > maxTotalDiss_ )
  {
    System::out() << myName_ << " says: TERMINATE because "
      << " totalDiss > maxTotalDiss" << endl;
    params.set ( XActions::TERMINATE, "sure" );
  }

  if ( ! isDispControl_ )
  {
    // don't accept solution and request switch
    // if the load changed sign

    if ( dispBC_->getLoadScale() * dispBC_->getLoadScale0() < 0. )
    {
      System::warn() << "Load changed sign.\n";

      params.set ( XActions::DISCARD, true );
    }

    // ... or if the displacement increment obtained with 
    // arclength is very large 

    else if ( dispBC_->getLastDispIncr( globdat ) / dispIncr0_ > 50. ) 
    {
      System::warn() << "Unacceptably large displacement increment.\n" ;
                     // << dispBC_->getLastDispIncr( globdat ) << " " 
                     // << dispIncr0_ << " (ignored).\n";  

      // params.set ( XActions::DISCARD, true );
    }

  }
}

//-----------------------------------------------------------------------
//   checkSwitch_
//-----------------------------------------------------------------------

void DispArclenModel::checkSwitch_

  ( const Properties&  params,
    const Properties&  globdat )

{
  bool   doCheck = false;

  // store number of iterations for adaptive stepping

  idx_t  iiter;

  params.get ( iiter, SolverInfo::ITER_COUNT );

  maxNIter_ = max ( maxNIter_, iiter );

  if ( isDispControl_ && ! tempStep_ )
  {
    // compute and store energy dissipation

    bool assem = globdat.find ( lastArcl_, XProps::FE_DISSIPATION );

    if ( ! assem )
    {
      lastArcl_ = getReleasedEnergy_( params, globdat );
    }

    // check whether solution is allowed (input dependent check)

    if ( isTmpDisp_ )
    {
      params.find ( doCheck, XActions::CHECK_BOUNDS );

      if ( doCheck && lastArcl_ > maxIncr_ )
      {
        System::warn() << "Energy increment too large. Solution discarded.\n"
          << "This check is optional. " << XActions::CHECK_BOUNDS << '\n';

        params.set ( XActions::DISCARD, true );
      }
    }

    // check for switch to arclength and if necessary initialize arcLength_

    if ( lastArcl_ > swtEner_ ||
        ( iiter > swtIter_ && lastArcl_ > minIncr_ ) ) 
    {
      params.set ( XActions::DO_SWITCH, "please" );

      arcLength_ = 0.1 * lastArcl_;
    }
  }
  else
  {
    lastArcl_ = arcLength_;
    if ( preferDisp_ )
    {
      // if positive load increment: switch back to disp control

      double lambda0 = dispBC_->getLoadScale0();
      double dlambda = dispBC_->getLoadScale () - lambda0;

      onceDown_ |= dlambda < 0.;

      // NB: not using onceDown flag,
      // rather only switch when small increments 

      if ( ( dlambda > 0. ) && ( arcLength_ <= swtEner_ ) )
      {
        out_ << "Loads increases, while preferDisp_, " 
             << "do switch if accepted...\n";

        params.set ( XActions::DO_SWITCH, "ifAccept" );

        onceDown_ = false;
      }
    }
  }
}

//-----------------------------------------------------------------------
//   getStepSize_
//-----------------------------------------------------------------------

void DispArclenModel::getStepSize_

  ( const Properties&  params,
    const Properties&  globdat )

{
  double incr = isDispControl_ ? 
                dispBC_->getDispIncr() :
                arcLength_ ;

  params.set ( XProps::STEP_SIZE, incr );
}

//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void DispArclenModel::setStepSize_

  ( const Properties&  params,
    const Properties&  globdat )

{
  double incr;
  params.get ( incr, XProps::STEP_SIZE );

  if ( isDispControl_ ) 
  {
    dispBC_->setDispIncr ( incr );
  }
  else
  {
    arcLength_ = incr; 
  }
}

//-----------------------------------------------------------------------
//   reduceStep_
//-----------------------------------------------------------------------

void DispArclenModel::reduceStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // Reduce step size if and return DONE when admissible

  if ( isDispControl_ )
  {
    double newIncr = dispBC_->reduceDispIncr ( reduction_ );

    if ( jem::numeric::abs( newIncr ) > minDispIncr_ )
    {
      out_ << "Displacement increment reduced to "
           << newIncr << "\n";

      params.set ( XActions::DONE, true );
    }
  }
  else
  {
    if ( arcLength_ > minIncr_ )
    {
      arcLength_ *= reduction_;

      out_ << "Path following parameter reduced to "
           << arcLength_ << "\n";

      params.set ( XActions::DONE, true );
    }
  }
}

//-----------------------------------------------------------------------
//   increaseStep_
//-----------------------------------------------------------------------

void DispArclenModel::increaseStep_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // Increase step size and if and return DONE when admissible
  // Adapt not from current step size but rather from the maximum
  // value that has been tried for this time step

  double maxTried;

  params.get ( maxTried, XProps::ADAPT_FROM );

  if ( isDispControl_ )
  {
    dispBC_->setDispIncr ( maxTried );

    double newIncr = dispBC_->reduceDispIncr ( 1./reduction_ );

    if ( newIncr / dispIncr0_ < 1. )
    {
      out_ << "Displacement increment increased from " << maxTried << 
              " to " << newIncr << "\n";

      params.set ( XActions::DONE, true );
    }
  }
  else
  {
    arcLength_ = maxTried;

    if ( arcLength_ < maxIncr_ )
    {
      arcLength_ /= reduction_;

      out_ << "Path following parameter increased to "
           << arcLength_ << "\n";

      params.set ( XActions::DONE, true );
    }
  }
}

//-----------------------------------------------------------------------
//   toArclControl_
//-----------------------------------------------------------------------

void DispArclenModel::toArclControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  bool accepted  = true;
  globdat.find   ( accepted, "var.accepted" );
  
  bool skipCheck = isTmpDisp_ &! accepted;

  if ( ! skipCheck )
  {
    if        ( arcLength_ < minIncr_ )
    {
      arcLength_ = minIncr_;
    }
    else if   ( arcLength_ > maxIncr_ )
    {
      arcLength_ = maxIncr_;
    }
  }

  isTmpDisp_ = isDispControl_ = false;

  out_ << "Switching to arclength control with dEnergy: " 
       << arcLength_ << endl;
}

//-----------------------------------------------------------------------
//   toDispControl_
//-----------------------------------------------------------------------

void DispArclenModel::toDispControl_

  ( const Properties&  params,
    const Properties&  globdat )

{
  isDispControl_ = true;

  bool temporary = false;

  params.find ( temporary, XActions::TEMPORARY );

  if ( temporary )
  {
    out_ << "Switching temporarily to disp control "
         << "fixing the displacement at " << dispBC_->getDispValue() << "\n";

    isTmpDisp_ = true;
  }
  else
  {
    dispBC_->setDispIncr ( dispIncr0_ );

    out_ << "Switching back to disp control "
         << "with disp incr: " << dispIncr0_ << "\n";
  }
}


//-----------------------------------------------------------------------
//   connect_
//-----------------------------------------------------------------------


void DispArclenModel::connect_ ()
{
  using jem::util::connect;

  connect ( dofs_->newSizeEvent,   this, & Self::dofsChanged_ );
  connect ( dofs_->newOrderEvent,  this, & Self::dofsChanged_ );

  dofs_->resetEvents ();
}


//-----------------------------------------------------------------------
//   dofsChanged_
//-----------------------------------------------------------------------


void DispArclenModel::dofsChanged_ ()
{
  updated_ = 0;
}


//-----------------------------------------------------------------------
//   consChanged_
//-----------------------------------------------------------------------


void DispArclenModel::consChanged_ ()
{
}

//-----------------------------------------------------------------------
//   getReleasedEnergy_
//-----------------------------------------------------------------------

// compute the energy released during one load step
// used to check for switch to use ArclenModule 

double DispArclenModel::getReleasedEnergy_

 ( const Properties&  params,
   const Properties&  globdat )

{
  if ( ! (updated_ & U_LOAD_) )
  {
    initLoad_ ( globdat );
  }

  // get master dof index with nonzero displacement
  // [currently DispControlModel deals only with one, 
  //  but this function is more general]

//   IdxVector idofs;
// 
//   idx_t nLoaded = dispBC_->getLoadedDofs ( idofs );

  idx_t nLoaded = idofs_.size();

  double jac11, lambda, dlambda, lambda0, dEnergy;
  Vector jac10, u, u0, du, tmp, fDiss;

  StateVector::get    ( u,  dofs_, globdat );
  StateVector::getOld ( u0, dofs_, globdat );

  Vector loadB ( select ( unitLoad_, idofs_ ) );
  Vector uB    ( select ( u,         idofs_ ) );
  Vector u0B   ( select ( u0,        idofs_ ) );

  du. resize ( nLoaded );
  tmp.resize ( nLoaded );

  // update displacement value in child

  if ( params.find ( lambda,  ArclenParams::LOAD_SCALE ) )
  {
    // when called in evalArcFunc 

    lambda0 = dispBC_->getLoadScale0();
    params.get ( jac10,   ArclenParams::JACOBIAN10 );

    dlambda = lambda - lambda0;
    jac11   = -0.5 * dot ( u0B, loadB );
    jac10   = 0.;

    select( jac10, idofs_ ) = 0.5 * lambda0 * loadB;

    // System::out() << "jac10 " << jac10[idofs_] << ", jac11 " << jac11 << endl;
    params.set ( ArclenParams::JACOBIAN11, jac11 );

    if ( ! isDispControl_ )
    {
      dispBC_->setDispValue ( uB[0] );
      dispBC_->setLoadScale ( lambda );
    }
  }
  else
  {
    // when called in checkSwitch

    lambda0 = dispBC_->getLoadScale0();
    dlambda = dispBC_->getLoadScale () - lambda0;
  }

  axpy ( du, uB, -1.0, u0B );

  energy0_  = 0.5 * dot ( u0B, loadB );
  energy0_ *= lambda0;

  tmp = lambda0 * du - dlambda * u0B;

  dEnergy   = 0.5 * dot ( tmp, loadB );

  return dEnergy;
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   declareDispArclenModel
//-----------------------------------------------------------------------


void declareDispArclenModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( DispArclenModel::TYPE_NAME,
                          & DispArclenModel::makeNew );
}
