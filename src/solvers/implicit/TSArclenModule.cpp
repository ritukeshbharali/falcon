
/*
 *  Copyright (C) 2019 DRG. All rights reserved.
 *
 *  This file is part of Jive, an object oriented toolkit for solving
 *  partial differential equations.
 *
 *  Commercial License Usage
 *
 *  This file may be used under the terms of a commercial license
 *  provided with the software, or under the terms contained in a written
 *  agreement between you and DRG. For more information contact DRG at
 *  http://www.dynaflow.com.
 *
 *  GNU Lesser General Public License Usage
 *
 *  Alternatively, this file may be used under the terms of the GNU
 *  Lesser General Public License version 2.1 or version 3 as published
 *  by the Free Software Foundation and appearing in the file
 *  LICENSE.LGPLv21 and LICENSE.LGPLv3 included in the packaging of this
 *  file. Please review the following information to ensure the GNU
 *  Lesser General Public License requirements will be met:
 *  https://www.gnu.org/licenses/lgpl.html and
 *  http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
 *  This file is part of Jive, an object oriented toolkit for
 *  solving partial differential equations.
 *
 *  Jive version: 3.0
 *  Date:         Fri 20 Dec 14:30:12 CET 2019
 */

/** @file TSArclenModule.cpp
 *  @brief Implements a time-step computing arc-length solver
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 10 January 2024
 *
 *  Usage: same as the Jive ArclenModule.
 */


#include <cmath>
#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/ArithmeticException.h>
#include <jem/base/IllegalInputException.h>
#include <jem/io/Writer.h>
#include <jem/util/None.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/util/utilities.h>
#include <jive/util/Globdat.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/ModelFactory.h>
#include <jive/app/ModuleFactory.h>
#include <jive/implict/Names.h>
#include <jive/implict/Globdat.h>
#include <jive/implict/ConHandler.h>
#include <jive/implict/SolverInfo.h>
#include <jive/implict/NonlinRunData.h>
#include <jive/implict/ArclenActions.h>

#include "TSArclenModule.h"
#include "util/XNames.h"


JEM_DEFINE_CLASS( jive::implict::TSArclenModule );


JIVE_BEGIN_PACKAGE( implict )


using jem::max;
using jem::isTiny;
using jem::newInstance;
using jem::Float;
using jem::System;
using jem::ArithmeticException;
using jem::IllegalInputException;
using jem::io::endl;
using jem::numeric::axpy;
using jive::model::Actions;
using jive::model::StateVector;


//=======================================================================
//   class TSArclenModule::RunData_
//=======================================================================


class TSArclenModule::RunData_ : public NonlinRunData
{
 public:

  typedef NonlinRunData   Super;
  typedef RunData_        Self;


  explicit inline         RunData_

    ( const String&         context );

  void                    initArcModel

    ( const String&         name,
      const Properties&     conf,
      const Properties&     props,
      const Properties&     globdat );

  bool                    commit

    ( const Properties&     globdat );

  void                    cancel

    ( const Properties&     globdat );


 public:

  Ref<Model>              arcModel;

  Matrix                  vbuf;

  double                  loadScale;
  double                  loadScale0;
  idx_t                   iterCount;

};


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


inline TSArclenModule::RunData_::RunData_ ( const String& ctx ) :

  Super ( ctx )

{
  loadScale  = 0.0;
  loadScale0 = 0.0;
  iterCount  = 0;
}


//-----------------------------------------------------------------------
//   initArcModel
//-----------------------------------------------------------------------


void TSArclenModule::RunData_::initArcModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jem::util::None;
  using jem::util::isNone;
  using jive::model::ModelFactory;
  using jive::model::ActionParams;

  if ( ! globdat.find( arcModel, Globdat::ARCLEN_MODEL ) )
  {
    Ref<Object>  obj;

    props.find ( obj, name );

    if ( ! obj || isNone( obj ) )
    {
      conf.set ( name, None::getInstance() );

      arcModel = model;
    }
    else
    {
      arcModel = ModelFactory::newInstance ( name,  conf,
                                             props, globdat );
    }
  }

  if ( arcModel != model )
  {
    Properties  params;

    arcModel->configure  ( props, globdat );
    arcModel->getConfig  ( conf,  globdat );
    arcModel->takeAction ( Actions::INIT, params, globdat );
  }
}


//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------


bool TSArclenModule::RunData_::commit ( const Properties& globdat )
{
  using jive::util::Globdat;
  using jive::model::ActionParams;

  Properties  params;

  bool        accept = true;

  params.set ( ArclenParams::IITER,          iterCount  );
  params.set ( ArclenParams::LOAD_SCALE,     loadScale  );
  params.set ( ArclenParams::OLD_LOAD_SCALE, loadScale0 );

  model->takeAction ( Actions::CHECK_COMMIT, params, globdat );

  if ( arcModel != model )
  {
    arcModel->takeAction ( Actions::CHECK_COMMIT, params, globdat );
  }

  params.find ( accept, ActionParams::ACCEPT );

  if ( accept )
  {
    params.erase ( ActionParams::ACCEPT );

    model->takeAction ( Actions::COMMIT, params, globdat );

    if ( arcModel != model )
    {
      arcModel->takeAction ( Actions::COMMIT, params, globdat );
    }

    Globdat    ::commitStep ( globdat );

    StateVector::updateOldOld ( dofs, globdat );
    StateVector::updateOld    ( dofs, globdat );

    loadScale0 = loadScale;
  }

  return accept;
}


//-----------------------------------------------------------------------
//   cancel
//-----------------------------------------------------------------------


void TSArclenModule::RunData_::cancel ( const Properties& globdat )
{
  using jive::util::Globdat;

  Properties  params;

  loadScale = fabs(loadScale0)/2.0;   // absolute value timestep/2
  iterCount = 0;

  Globdat    ::restoreStep ( globdat );
  StateVector::restoreNew  ( dofs, globdat );

  params.set ( ArclenParams::IITER,          iterCount  );
  params.set ( ArclenParams::LOAD_SCALE,     loadScale  );
  params.set ( ArclenParams::OLD_LOAD_SCALE, loadScale0 );

  model->takeAction ( Actions::CANCEL, params, globdat );

  if ( arcModel != model )
  {
    arcModel->takeAction ( Actions::CANCEL, params, globdat );
  }
}


//=======================================================================
//   class TSArclenModule::Work_
//=======================================================================


class TSArclenModule::Work_ : public ConHandler
{
 public:

  typedef ArclenParams    Params;


                          Work_

    ( TSArclenModule&         module,
      const Properties&     globdat );

  inline                 ~Work_           ();

  void                    updateState     ();

  void                    updateRscale

    ( const Properties&     globdat );

  void                    updateResidual  ();

  inline void             updateModel

    ( const Properties&     globdat );

  void                    updateUnitLoad

    ( const Properties&     globdat );

  void                    updateArcFunc

    ( const Properties&     globdat );

  bool                    checkConverged

    ( double                tol,
      const Properties&     globdat );

  void                    reportProgress  ();

  void                    updateStepSize

    ( const Properties&     globdat );


 public:

  RunData_&               rundat;

  Vector                  u;
  Vector                  d;
  Vector                  q;
  Vector                  r;
  Vector                  fint;
  Vector                  fext;
  Vector                  load;

  Vector                  jac10;
  double                  jac11;
  double                  arcFunc;
  double                  loadScale;

  double                  rnorm;
  double                  rscale;

  idx_t                   iiter;

  Properties              params;

 private:

  String                  myName_;

};


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


TSArclenModule::Work_::Work_

  ( TSArclenModule&      mod,
    const Properties&  globdat ) :

    rundat  ( *mod.rundat_ ),
    myName_ (  mod.getName() )

{
  idx_t  dofCount;
  idx_t  j;


  jac11     = 0.0;
  arcFunc   = 0.0;
  loadScale = rundat.loadScale;
  rnorm     = 0.0;
  rscale    = 0.0;
  iiter     = 0;

  rundat.updateConstraints ( globdat );
  rundat.dofs->resetEvents ();

  rundat.frozen = true;
  dofCount      = rundat.dofs->dofCount ();

  StateVector::get ( u, rundat.dofs, globdat );

  j = 0;

  rundat.vbuf.resize ( dofCount, 7 );

  d    .ref ( rundat.vbuf[j++] );
  q    .ref ( rundat.vbuf[j++] );
  r    .ref ( rundat.vbuf[j++] );
  fint .ref ( rundat.vbuf[j++] );
  fext .ref ( rundat.vbuf[j++] );
  load .ref ( rundat.vbuf[j++] );
  jac10.ref ( rundat.vbuf[j++] );

  params.set ( Params::IITER,          iiter );
  params.set ( Params::UNIT_LOAD,      load );
  params.set ( Params::LOAD_SCALE,     loadScale );
  params.set ( Params::OLD_LOAD_SCALE, rundat.loadScale0 );
  params.set ( Params::LOAD_RESPONSE,  q );
  params.set ( Params::JACOBIAN10,     jac10 );

  saveConstraints ( *rundat.cons );

  if ( ! (mod.options_ & DELTA_CONS) )
  {
    adjustConstraints ( *rundat.cons, u );
  }
}


inline TSArclenModule::Work_::~Work_ ()
{
  rundat.frozen    = false;
  rundat.loadScale = loadScale;
  rundat.iterCount = iiter;
}


//-----------------------------------------------------------------------
//   updateState
//-----------------------------------------------------------------------


void TSArclenModule::Work_::updateState ()
{
  double  alpha, b, c;


  b = arcFunc + rundat.vspace->product ( jac10, d );
  c = jac11   + rundat.vspace->product ( jac10, q );

  if ( Float::isNaN( b + c ) )
  {
    throw ArithmeticException (
      rundat.context,
      "error evaluating the arc-length function: NaN"
    );
  }

  if ( isTiny( c ) )
  {
    throw ArithmeticException (
      rundat.context,
      "singular load increment"
    );
  }

  alpha = b / c;

  axpy ( u,    1.0, d );
  axpy ( u, -alpha, q );

  loadScale -= alpha;
  iiter++;

  params.set ( Params::IITER,      iiter     );
  params.set ( Params::LOAD_SCALE, loadScale );
}


//-----------------------------------------------------------------------
//   updateRscale
//-----------------------------------------------------------------------


void TSArclenModule::Work_::updateRscale

  ( const Properties&  globdat )

{
  using jive::model::ActionParams;

  Properties  params;
  double      rtmp;

  r      = fext - fint;

  rtmp   = std::sqrt( Float::EPSILON ) *

    max ( rundat.vspace->norm2( fint ),
          rundat.vspace->norm2( fext ) );

  rtmp   = max ( rtmp, rundat.vspace->norm2( r ) );
  rscale = max ( rtmp, rscale );

  params.set ( ActionParams::RESIDUAL,   r );
  params.set ( ActionParams::RES_SCALE,  rscale );
  params.set ( ActionParams::INT_VECTOR, fint );
  params.set ( ActionParams::EXT_VECTOR, fext );

  rundat.model->takeAction ( Actions::GET_RES_SCALE,
                             params, globdat );

  params.get ( rscale, ActionParams::RES_SCALE );
}


//-----------------------------------------------------------------------
//   updateResidual
//-----------------------------------------------------------------------


void TSArclenModule::Work_::updateResidual ()
{
  r = fext-fint;

  //axpy           ( r, 1.0, load );
  evalMasterDofs ( r, *rundat.cons );

  rnorm = rundat.vspace->norm2 ( r );

  if ( rscale > 0.0 )
  {
    rnorm /= rscale;
  }
}


//-----------------------------------------------------------------------
//   updateModel
//-----------------------------------------------------------------------


inline void TSArclenModule::Work_::updateModel

  ( const Properties&  globdat )

{
  rundat.model->takeAction ( Actions::UPDATE, params, globdat );

  if ( rundat.arcModel != rundat.model )
  {
    rundat.arcModel->takeAction ( Actions::UPDATE, params, globdat );
  }
}


//-----------------------------------------------------------------------
//   updateUnitLoad
//-----------------------------------------------------------------------


void TSArclenModule::Work_::updateUnitLoad

  ( const Properties&  globdat )

{
  double  tmp;
  bool    result;


  load   = 0.0;
  result = rundat.arcModel ->

    takeAction ( ArclenActions::GET_UNIT_LOAD, params, globdat );

  if ( ! result )
  {
    String  name = rundat.arcModel->getName ();

    throw IllegalInputException (
      rundat.context,
      String::format (
        "no load returned by the model `%s\' "
        "(ignored action `%s\')",
        name,
        ArclenActions::GET_UNIT_LOAD
      )
    );
  }

  tmp = rundat.vspace->norm2 ( load );

  if ( Float::isNaN( tmp ) )
  {
     throw ArithmeticException (
      rundat.context,
      "load vector contains invalid data: NaN"
     );
  }

  if ( isTiny( tmp ) )
  {
    print ( System::warn(), rundat.context,
            " : zero load vector.\n" );
  }
}


//-----------------------------------------------------------------------
//   updateArcFunc
//-----------------------------------------------------------------------


void TSArclenModule::Work_::updateArcFunc

  ( const Properties&  globdat )

{
  using jive::util::evalMasterDofs;

  bool  result;


  jac10  = 0.0;
  jac11  = 0.0;
  result = rundat.arcModel ->

    takeAction ( ArclenActions::GET_ARC_FUNC, params, globdat );

  if ( ! result )
  {
    String  name = rundat.arcModel->getName ();

    throw jem::Exception (
      rundat.context,
      String::format (
        "no arc-length function returned by the model `%s\' "
        "(ignored action `%s\')",
        name,
        ArclenActions::GET_ARC_FUNC
      )
    );
  }

  params.get     ( arcFunc, ArclenParams::ARC_FUNC   );
  params.find    ( jac11,   ArclenParams::JACOBIAN11 );

  evalMasterDofs ( jac10, *rundat.cons );
}


//-----------------------------------------------------------------------
//   checkConverged
//-----------------------------------------------------------------------


bool TSArclenModule::Work_::checkConverged

  ( double             tol,
    const Properties&  globdat )

{
  using jive::model::ActionParams;

  bool  conv = (rnorm <= tol && std::fabs( arcFunc ) <= tol);

  params.set ( ActionParams::RESIDUAL,  r );
  params.set ( ActionParams::RES_SCALE, rscale );
  params.set ( ActionParams::CONVERGED, conv );

  rundat.model->takeAction ( Actions::CHECK_CONVERGED,
                             params, globdat );

  if ( rundat.arcModel != rundat.model )
  {
    rundat.arcModel->takeAction ( Actions::CHECK_CONVERGED,
                                  params, globdat );
  }

  params.get ( conv, ActionParams::CONVERGED );

  return conv;
}


//-----------------------------------------------------------------------
//   reportProgress
//-----------------------------------------------------------------------


void TSArclenModule::Work_::reportProgress ()
{
  NumberFormat&  nf = rundat.nformat;

  double   loadIncr = loadScale - rundat.loadScale0;


  print ( System::info( myName_ ),
          rundat.context, ", iteration ", iiter, endl );

  print ( System::info( myName_ ),
          "  load scale      = ", nf.print( loadScale ), endl );
  print ( System::info( myName_ ),
          "  load increment  = ", nf.print( loadIncr ),  endl );
  print ( System::info( myName_ ),
          "  arc-length func = ", nf.print( arcFunc ),   endl );
  print ( System::info( myName_ ),
          "  scaled residual = ", nf.print( rnorm ),     endl );

  System::info( myName_ ).flush ();
}


//-----------------------------------------------------------------------
//   updateStepSize
//-----------------------------------------------------------------------


void TSArclenModule::Work_::updateStepSize

  ( const Properties&  globdat )

{
  System::out() << "TimeStepArclen:: tstep = " << loadScale << "\n";

  params.set ( XProps::STEP_SIZE, loadScale );

  bool    result;

  result = rundat.model ->

    takeAction ( XActions::SET_STEP_SIZE, params, globdat );

  if ( ! result )
  {
    String  name = rundat.model->getName ();

    throw IllegalInputException (
      rundat.context,
      String::format (
        "model `%s\' did not respond to "
        "action `%s\'",
        name,
        XActions::SET_STEP_SIZE
      )
    );
  }

}


//=======================================================================
//   class TSArclenModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  TSArclenModule::TYPE_NAME  = "TimeStepArclen";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


TSArclenModule::TSArclenModule ( const String& name ) :

  Super ( name )

{
  maxIter_   = 20;
  options_   = 0;
  precision_ = 1.0e-3;
  loadScale_ = 0.0;
}


TSArclenModule::~TSArclenModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status TSArclenModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jive::util::joinNames;

  rundat_ = nullptr;

  Ref<RunData_>  newdat = newInstance<RunData_> ( getContext() );
  String         name   = joinNames ( myName_, PropNames::SOLVER );

  newdat->loadScale     = loadScale_;
  newdat->loadScale0    = loadScale_;

  newdat->init         ( globdat );
  newdat->initSolver   ( name, precision_, conf, props, globdat );

  name = joinNames     ( myName_, PropNames::MODEL );

  newdat->initArcModel ( name, conf, props, globdat );

  // Everything OK, so commit the changes

  rundat_.swap ( newdat );

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void TSArclenModule::shutdown ( const Properties& globdat )
{
  rundat_ = nullptr;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void TSArclenModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::maxOf;

  if ( props.contains( myName_ ) )
  {
    Properties  myProps = props.findProps ( myName_ );

    double      scale;
    bool        option;


    if ( myProps.find( option, PropNames::DELTA_CONS ) )
    {
      if ( option )
      {
        options_ |=  DELTA_CONS;
      }
      else
      {
        options_ &= ~DELTA_CONS;
      }
    }

    myProps.find ( maxIter_,   PropNames::MAX_ITER,
                   0,          maxOf( maxIter_ ) );
    myProps.find ( precision_, PropNames::PRECISION,
                   0.0,        1.0e20 );

    if ( myProps.find( scale, PropNames::LOAD_SCALE ) )
    {
      setLoadScale ( scale );
    }
  }

  if ( rundat_ )
  {
    RunData_&  d = * rundat_;

    d.solver->configure ( props );

    if ( d.arcModel != d.model )
    {
      d.arcModel->configure ( props, globdat );
    }
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void TSArclenModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  myConf.set ( PropNames::DELTA_CONS,
               (bool) (options_ & DELTA_CONS) );

  myConf.set ( PropNames::MAX_ITER,   maxIter_ );
  myConf.set ( PropNames::PRECISION,  precision_ );
  myConf.set ( PropNames::LOAD_SCALE, getLoadScale() );

  if ( rundat_ )
  {
    RunData_&  d = * rundat_;

    d.solver->getConfig ( conf );

    if ( d.arcModel != d.model )
    {
      d.arcModel->getConfig ( conf, globdat );
    }
  }
}


//-----------------------------------------------------------------------
//   advance
//-----------------------------------------------------------------------


void TSArclenModule::advance ( const Properties& globdat )
{
  if ( ! rundat_ )
  {
    notAliveError ( JEM_FUNC );
  }

  rundat_->advance ( globdat );
}


//-----------------------------------------------------------------------
//   solve
//-----------------------------------------------------------------------


void TSArclenModule::solve

  ( const Properties&  info,
    const Properties&  globdat )

{
  if ( ! rundat_ )
  {
    notAliveError ( JEM_FUNC );
  }

  RunData_&  d = * rundat_;

  Work_      w ( *this, globdat );

  double     loadIncr;
  bool       result;


  print ( System::info( myName_ ),
          "Starting the time-step based arc-length solver `",
          myName_, "\' ...\n" );

  globdat.set("iiter", w.iiter);

  d.getExtVector    ( w.fext, globdat );
  w.updateStepSize  ( globdat );
  w.updateUnitLoad  ( globdat );
  d.updateMatrix    ( w.fint, globdat );
  w.updateRscale    ( globdat );
  d.solver->solve   ( w.d, w.r );
  w.zeroConstraints ( *d.cons );
  d.solver->solve   ( w.q, w.load );
  w.updateArcFunc   ( globdat );

  result = solve_   ( w, globdat );

  w.restoreConstraints ( *d.cons );

  loadIncr = w.loadScale - d.loadScale0;

  info.set ( SolverInfo::CONVERGED,  result      );
  info.set ( SolverInfo::ITER_COUNT, w.iiter     );
  info.set ( SolverInfo::RESIDUAL,   w.rnorm     );
  info.set ( SolverInfo::LOAD_INCR,  loadIncr    );
  info.set ( SolverInfo::LOAD_SCALE, w.loadScale );

  if ( Globdat::hasVariable( myName_, globdat ) )
  {
    Properties  vars = Globdat::getVariables ( myName_, globdat );

    vars.set ( SolverInfo::CONVERGED,  result      );
    vars.set ( SolverInfo::ITER_COUNT, w.iiter     );
    vars.set ( SolverInfo::RESIDUAL,   w.rnorm     );
    vars.set ( SolverInfo::LOAD_INCR,  loadIncr    );
    vars.set ( SolverInfo::LOAD_SCALE, w.loadScale );
  }

  if ( ! result )
  {
    throw jem::Exception (
      getContext (),
      String::format (
        "no convergence achieved in %d iterations; "
        "final residual: %e",
        w.iiter,
        w.rnorm
      )
    );
  }

  print ( System::info( myName_ ),
          "The generalized arc-length solver converged in ",
          w.iiter, " iterations\n\n" );
}


//-----------------------------------------------------------------------
//   cancel
//-----------------------------------------------------------------------


void TSArclenModule::cancel ( const Properties& globdat )
{
  if ( ! rundat_ )
  {
    notAliveError ( JEM_FUNC );
  }

  rundat_->cancel ( globdat );
}


//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------


bool TSArclenModule::commit ( const Properties& globdat )
{
  if ( ! rundat_ )
  {
    notAliveError ( JEM_FUNC );
  }

  return rundat_->commit ( globdat );
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void TSArclenModule::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps > 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double TSArclenModule::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void TSArclenModule::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void TSArclenModule::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


TSArclenModule::Options TSArclenModule::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setLoadScale
//-----------------------------------------------------------------------


void TSArclenModule::setLoadScale ( double scale )
{
  if ( ! rundat_ )
  {
    loadScale_ = scale;
  }
  else
  {
    rundat_->loadScale  = scale;
    rundat_->loadScale0 = scale;
  }
}


//-----------------------------------------------------------------------
//   getLoadScale
//-----------------------------------------------------------------------


double TSArclenModule::getLoadScale () const
{
  if ( ! rundat_ )
  {
    return loadScale_;
  }
  else
  {
    return rundat_->loadScale;
  }
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Module> TSArclenModule::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  return newInstance<Self> ( name );
}


//-----------------------------------------------------------------------
//   declare
//-----------------------------------------------------------------------


void TSArclenModule::declare ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( TYPE_NAME,  & makeNew );
  ModuleFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   solve_
//-----------------------------------------------------------------------


bool TSArclenModule::solve_

  ( Work_&             w,
    const Properties&  globdat )

{
  RunData_&  d = w.rundat;


  w.updateState    ();
  w.updateModel    ( globdat );
  w.updateStepSize ( globdat );
  w.updateUnitLoad ( globdat );
  w.updateArcFunc  ( globdat );
  d.updateMatrix   ( w.fint, globdat );
  w.updateRscale   ( globdat );

  print ( System::info( myName_ ), d.context,
          " : residual scale factor = ",
          d.nformat.print( w.rscale ), endl );

  w.updateResidual ();

  do
  {
    w.reportProgress ();

    if ( w.checkConverged( precision_, globdat ) )
    {
      return true;
    }

    if ( w.iiter > maxIter_ )
    {
      return false;
    }

    if ( w.rnorm > 1.0e4 || std::fabs( w.arcFunc ) > 1.0e4 )
    {
      return false;
    }

    d.solver->solve ( w.d, w.r    );
    d.solver->solve ( w.q, w.load );

    w.updateState    ();
    w.updateModel    ( globdat );
    w.updateStepSize ( globdat );
    w.updateUnitLoad ( globdat );
    w.updateArcFunc  ( globdat );
    d.updateMatrix   ( w.fint, globdat );
    w.updateResidual ();
  }
  while ( true );

  return false;
}


JIVE_END_PACKAGE( implict )
