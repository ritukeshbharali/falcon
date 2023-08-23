
#include <cmath>
#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/System.h>
#include <jem/base/Exception.h>
#include <jem/base/ClassTemplate.h>
#include <jem/base/array/operators.h>
#include <jem/io/Writer.h>
#include <jem/util/Event.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/util/utilities.h>
#include <jive/util/Globdat.h>
#include <jive/algebra/VectorSpace.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/app/ModuleFactory.h>
#include <jive/implict/Names.h>
#include <jive/implict/ConHandler.h>
#include <jive/implict/SolverInfo.h>
#include <jive/implict/SolverBounds.h>
#include <jive/implict/NonlinRunData.h>

#include "BFGSModule.h"
#include "util/XNames.h"


JEM_DEFINE_CLASS( jive::implict::BFGSModule );


JIVE_BEGIN_PACKAGE( implict )


using jem::min;
using jem::max;
using jem::isTiny;
using jem::newInstance;
using jem::Float;
using jem::System;
using jem::Exception;
using jem::io::endl;
using jem::numeric::axpy;
using jive::model::Actions;
using jive::model::StateVector;


/*
 * The BFGSModule implements both the traditional Newton-Raphson
 * method and the quasi-newton BFGS method; this is determined by
 * the REFORM_ITER option. When the latter is set to a positive
 * value then the BFGS method is applied.
 *
 * The BFGS method has been implemented by Yuli Huang (Arup) and
 * Jian-Ying Wu (SCUT). The method is described in:
 *
 *   J.Y. Wu, Y. Huang and V.P. Nguyen. "On the BFGS
 *   monolithic algorithm for the unified phase-field damage
 *   theory". Computer Methods in Applied Mechanics and
 *   Engineering, 112704, 2019.
 *
 *   J.Y. Wu and Y Huang. "Comprehensive ABAQUS implementation of
 *   phase-field damage models for fracture in solids".
 *   Theoretical and Applied Fracture Mechanics, 2019, in press.
 *
 * The BFGSModule can solve generic non-linear systems of equations
 * and non-linear complementary problems with box constraints. The
 * latter type of problems arise in minimisation problems in which
 * the unknowns (DOFs) must remain within specific bounds.
 *
 * See also:
 *
 *   Steven J Benson and Todd S Munson. "Flexible complementary
 *   solvers for large-scale applications". Optimization Methods
 *   and Software, Vol 21, No 1, February 2006, pp 155-168.
 *
 *   Francisco Facchinei, Joaquim Judice, and Joao Soares. "An
 *   active-set Newton algorithm for large-scale non-linear
 *   programs with box constraints". SIAM J Optim., Vol 8,
 *   No 1, February 1998, pp 158-186.
 */

//=======================================================================
//   class BFGSModule::RunData_
//=======================================================================


class BFGSModule::RunData_ : public NonlinRunData
{
 public:

  typedef NonlinRunData   Super;
  typedef RunData_        Self;


  explicit inline         RunData_

    ( const String&         context );


 public:

  Matrix                  vbuf;
  IdxMatrix               ibuf;
  Matrix                  pvecs;
  Matrix                  qvecs;

  bool                    validMatrix;


 protected:

  virtual void            dofsChanged_  () override;

};


//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


inline BFGSModule::RunData_::RunData_ ( const String& ctx ) :

  Super ( ctx )

{
  validMatrix = false;
}


//-----------------------------------------------------------------------
//   dofsChanged_
//-----------------------------------------------------------------------


void BFGSModule::RunData_::dofsChanged_ ()
{
  if ( frozen )
  {
    Super::dofsChanged_ ();
  }
  else
  {
    validMatrix = false;
  }
}


//=======================================================================
//   class BFGSModule::Work_
//=======================================================================


class BFGSModule::Work_ : public ConHandler
{
 public:

                          Work_

    ( BFGSModule&         module,
      const Properties&     globdat );

  inline                 ~Work_             ();

  void                    updateBounds

    ( const Properties&     globdat );

  void                    applyBounds       ();

  inline bool             hasBounds         () const noexcept;

  void                    updateRscale

    ( const Properties&     globdat );

  void                    updateResidual    ();

  void                    updateMatrix

    ( const Properties&     globdat );

  void                    calcIncrement

    ( double                maxIncr );

  void                    doTotalUpdate

    ( double                scale,
      const Properties&     globdat );

  bool                    checkConverged

    ( double                tol,
      const Properties&     globdat );

  inline void             reportProgress    ();


 public:

  RunData_&               rundat;

  Vector                  u;
  Vector                  u0;
  Vector                  du;
  Vector                  fext;
  Vector                  fint;
  Vector                  r;
  Vector                  r0;
  Vector                  r2;
  Vector                  dr;

  Matrix                  pvecs;
  Matrix                  qvecs;

  IdxVector               slaveDofs;

  idx_t                   iiter;
  idx_t                   jiter;
  double                  rscale;
  double                  uscale;
  double                  dnorm;
  double                  rnorm;
  double                  enorm;


 private:

  String                  myName_;
  Vector                  lowerBound_;
  Vector                  upperBound_;
  IdxVector               fixedDofs_;
  bool                    bounded_;

};


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


BFGSModule::Work_::Work_

  ( BFGSModule&      mod,
    const Properties&  globdat ) :

    rundat  ( *mod.rundat_ ),
    myName_ (  mod.getName() )

{
  idx_t  dofCount;
  idx_t  j;


  rundat.updateConstraints ( globdat );
  rundat.dofs->resetEvents ();

  iiter              = 0;
  jiter              = 0;
  rscale             = 0.0;
  uscale             = 1.0;
  dnorm              = 0.0;
  rnorm              = 1.0;
  enorm              = 0.0;
  bounded_           = mod.bounds_->hasBounds ();
  rundat.validMatrix = false;
  rundat.frozen      = true;
  dofCount           = rundat.dofs->dofCount ();

  StateVector::get ( u, rundat.dofs, globdat );

  j = 0;

  rundat.vbuf .resize ( dofCount, 10 );
  rundat.ibuf .resize ( dofCount,  1 );
  rundat.pvecs.resize ( dofCount, mod.reformIter_ );
  rundat.qvecs.resize ( dofCount, mod.reformIter_ );

  u0         .ref ( rundat.vbuf[j++] );
  du         .ref ( rundat.vbuf[j++] );
  fext       .ref ( rundat.vbuf[j++] );
  fint       .ref ( rundat.vbuf[j++] );
  r          .ref ( rundat.vbuf[j++] );
  r0         .ref ( rundat.vbuf[j++] );
  r2         .ref ( rundat.vbuf[j++] );
  dr         .ref ( rundat.vbuf[j++] );
  lowerBound_.ref ( rundat.vbuf[j++] );
  upperBound_.ref ( rundat.vbuf[j++] );

  pvecs      .ref ( rundat.pvecs );
  qvecs      .ref ( rundat.qvecs );

  j = 0;

  fixedDofs_ .ref ( rundat.ibuf[j++] );

  mod.bounds_->getBounds ( lowerBound_,
                           upperBound_,
                           *rundat.dofs );

  saveConstraints ( *rundat.cons );

  if ( ! (mod.options_ & DELTA_CONS) )
  {
    adjustConstraints ( *rundat.cons, u );
  }
}


inline BFGSModule::Work_::~Work_ ()
{
  rundat.frozen = false;
}


//-----------------------------------------------------------------------
//   updateBounds
//-----------------------------------------------------------------------


void BFGSModule::Work_::updateBounds

  ( const Properties&  globdat )

{
  bool  result = rundat.getBounds ( lowerBound_,
                                    upperBound_,
                                    globdat );

  bounded_ = bounded_ || result;
}


//-----------------------------------------------------------------------
//   applyBounds
//-----------------------------------------------------------------------


void BFGSModule::Work_::applyBounds ()
{
  if ( bounded_ )
  {
    const idx_t  dofCount = u.size ();

    if ( u          .isContiguous() &&
         lowerBound_.isContiguous() &&
         upperBound_.isContiguous() )
    {
      const double* JEM_RESTRICT  lb = lowerBound_.addr ();
      const double* JEM_RESTRICT  ub = upperBound_.addr ();
      double*       JEM_RESTRICT  x  = u          .addr ();

      for ( idx_t i = 0; i < dofCount; i++ )
      {
        if      ( x[i] < lb[i] )
        {
          x[i] = lb[i];
        }
        else if ( x[i] > ub[i] )
        {
          x[i] = ub[i];
        }
      }
    }
    else
    {
      for ( idx_t i = 0; i < dofCount; i++ )
      {
        if      ( u[i] < lowerBound_[i] )
        {
          u[i] = lowerBound_[i];
        }
        else if ( u[i] > upperBound_[i] )
        {
          u[i] = upperBound_[i];
        }
      }
    }
  }
}


//-----------------------------------------------------------------------
//   hasBounds()
//-----------------------------------------------------------------------


inline bool BFGSModule::Work_::hasBounds () const noexcept
{
  return bounded_;
}


//-----------------------------------------------------------------------
//   updateRscale
//-----------------------------------------------------------------------


void BFGSModule::Work_::updateRscale

  ( const Properties&  globdat )

{
  using jive::model::ActionParams;

  Properties  params;
  double      rtmp;

  axpy ( r, fext, -1.0, fint );

  rtmp   = std::sqrt( Float::EPSILON ) *

    max ( rundat.vspace->norm2( fext ),
          rundat.vspace->norm2( fint ) );

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


void BFGSModule::Work_::updateResidual ()
{
  using std::fabs;

  axpy           ( r, fext, -1.0, fint );
  evalMasterDofs ( r, *rundat.cons );

  if ( ! bounded_ )
  {
    rnorm = rundat.vspace->norm2   ( r );
    enorm = rundat.vspace->product ( r, du );
  }
  else
  {
    const double  epsilon  = std::sqrt ( Float::EPSILON );
    const idx_t   dofCount = u.size    ();

    for ( idx_t i = 0; i < dofCount; i++ )
    {
      double  lb =  lowerBound_[i];
      double  ub =  upperBound_[i];
      double  x  =  u[i];
      double  f  = -r[i];

      lb += epsilon * fabs ( lb );
      ub -= epsilon * fabs ( ub );

      // Note that the condition ((x <= lb) && (x >= ub)) may be true
      // if the lower bound equals the upper bound.

      if ( x <= lb )
      {
        f = min ( 0.0, f );
      }

      if ( x >= ub )
      {
        f = max ( 0.0, f );
      }

      r2[i] = f;
    }

    rnorm = rundat.vspace->norm2   ( r2 );
    enorm = rundat.vspace->product ( r2, du );
  }

  if ( rscale > 0.0 )
  {
    rnorm /= rscale;
  }
}


//-----------------------------------------------------------------------
//   updateMatrix
//-----------------------------------------------------------------------


void BFGSModule::Work_::updateMatrix ( const Properties& globdat )
{
  if ( ! rundat.validMatrix || (jiter >= pvecs.size(1)) )
  {
    rundat.updateMatrix ( fint, globdat );

    this ->jiter       = 0;
    rundat.validMatrix = true;

    print ( System::debug( myName_ ), rundat.context,
            " : global (tangent) matrix updated\n" );
  }
  else
  {
    rundat.getIntVector ( fint, globdat );
  }
}


//-----------------------------------------------------------------------
//   calcIncrement
//-----------------------------------------------------------------------


void BFGSModule::Work_::calcIncrement ( double maxIncr )
{
  using std::fabs;

  idx_t   fixedCount = 0;
  double  dnorm0     = dnorm;


  if ( bounded_ )
  {
    const double  epsilon  = std::sqrt ( Float::EPSILON );
    const idx_t   dofCount = u.size    ();

    for ( idx_t i = 0; i < dofCount; i++ )
    {
      double  lb =  lowerBound_[i];
      double  ub =  upperBound_[i];
      double  x  =  u[i];
      double  f  = -r[i];

      lb += epsilon * fabs ( lb );
      ub -= epsilon * fabs ( ub );

      if ( ((x <= lb) && (f > 0.0)) ||
           ((x >= ub) && (f < 0.0)) )
      {
        if ( ! rundat.cons->isSlaveDof( i ) )
        {
          fixedDofs_[fixedCount++] = i;

          rundat.cons->addConstraint ( i );
        }
      }
    }
  }

  r0 = r;
  u0 = u;

  for ( idx_t j = jiter - 1_idx; j >= 0; j-- )
  {
    axpy ( r, -rundat.vspace->product( qvecs[j], r ), pvecs[j] );
  }

  if ( ! bounded_ )
  {
    rundat.solver->solve ( du, r );
  }
  else
  {
    try
    {
      rundat.solver->solve ( du, r );
    }
    catch ( const Exception& )
    {
      // Apply a steepest descent step.

      print ( System::info( myName_ ), rundat.context,
              " : singular matrix; falling back "
              "to steepest descent\n" );

      du = r;
    }
  }

  for ( idx_t j = 0 ; j < jiter; j++ )
  {
    axpy ( du, -rundat.vspace->product( pvecs[j], du ), qvecs[j] );
  }

  for ( idx_t i = 0 ; i < fixedCount; i++ )
  {
    rundat.cons->eraseConstraint ( fixedDofs_[i] );
  }

  dnorm = rundat.vspace->norm2 ( du );

  iiter++;

  if ( bounded_ || (iiter <= 2) || isTiny( dnorm0 ) )
  {
    return;
  }

  if ( dnorm > maxIncr * dnorm0 )
  {
    double  scale = maxIncr * dnorm0 / dnorm;

    print ( System::info( myName_ ), rundat.context,
            " : scaling solution increment with factor ",
            rundat.nformat.print( scale ),
            endl );

    du    *= scale;
    dnorm *= scale;
  }
}


//-----------------------------------------------------------------------
//   doTotalUpdate
//-----------------------------------------------------------------------


void BFGSModule::Work_::doTotalUpdate

  ( double             scale,
    const Properties&  globdat )

{
  uscale = scale;

  axpy                ( u, u0, scale, du );
  applyBounds         ();
  rundat.updateModel  (       globdat );
  rundat.getIntVector ( fint, globdat );
  updateResidual      ();
}


//-----------------------------------------------------------------------
//   checkConverged
//-----------------------------------------------------------------------


bool BFGSModule::Work_::checkConverged

  ( double             tol,
    const Properties&  globdat )

{
  using jive::model::ActionParams;

  Properties  params;
  bool        conv;


  conv = (rnorm <= tol);

  params.set ( ActionParams::RESIDUAL,  r );
  params.set ( ActionParams::RES_SCALE, rscale );
  params.set ( ActionParams::CONVERGED, conv );

  rundat.model->takeAction ( Actions::CHECK_CONVERGED,
                             params, globdat );

  params.find ( conv, ActionParams::CONVERGED );

  return conv;
}


//-----------------------------------------------------------------------
//   reportProgress
//-----------------------------------------------------------------------


inline void BFGSModule::Work_::reportProgress ()
{
  print ( System::info( myName_ ), rundat.context,
          " : iter = "          , iiter,
          ", scaled residual = ", rundat.nformat.print( rnorm ),
          endl );

  System::info( myName_ ).flush ();
}


//=======================================================================
//   class BFGSModule
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  BFGSModule::TYPE_NAME = "BFGS";


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


BFGSModule::BFGSModule ( const String& name ) :

  Super ( name )

{
  maxIter_    = 20;
  reformIter_ = 0;
  options_    = 0;
  tiny_       = jem::Limits<double>::TINY_VALUE;
  precision_  = 1.0e-3;
  lsearchTol_ = 0.9;
  maxIncr_    = 10.0;
  maxResIncr_ = 0.99;
  bounds_     = jem::newInstance<SolverBounds> ( getContext() );
}


BFGSModule::~BFGSModule ()
{}


//-----------------------------------------------------------------------
//   init
//-----------------------------------------------------------------------


Module::Status BFGSModule::init

  ( const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat )

{
  using jive::util::joinNames;

  rundat_ = nullptr;

  Ref<RunData_>  newdat = newInstance<RunData_> ( getContext() );
  String         name   = joinNames ( myName_, PropNames::SOLVER );

  newdat->init       ( globdat );
  newdat->initSolver ( name, precision_, conf, props, globdat );

  // Everything OK, so commit the changes

  rundat_.swap ( newdat );

  return OK;
}


//-----------------------------------------------------------------------
//   shutdown
//-----------------------------------------------------------------------


void BFGSModule::shutdown ( const Properties& globdat )
{
  rundat_ = nullptr;
}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void BFGSModule::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  using jem::maxOf;

  if ( props.contains( myName_ ) )
  {
    Properties  myProps = props.findProps ( myName_ );

    bool        option;

    myProps.find ( maxIter_,    PropNames::MAX_ITER,
                   0,           maxOf( maxIter_ ) );
    myProps.find ( reformIter_, XProps::REFORM_ITER,
                   0,           maxIter_ );
    myProps.find ( tiny_,       PropNames::TINY,
                   0.0,         1.0e20 );
    myProps.find ( precision_,  PropNames::PRECISION,
                   0.0,         1.0e20 );
    myProps.find ( lsearchTol_, XProps::LINE_SEARCH_TOL,
                   0.0,         1.0e20 );
    myProps.find ( maxIncr_,    PropNames::MAX_INCR,
                   0.0,         1.0e20 );
    myProps.find ( maxResIncr_, XProps::MAX_RES_INCR,
                   0.0,         1.0e20 );

    if ( myProps.find( option, PropNames::LINE_SEARCH ) )
    {
      if ( option )
      {
        options_ |=  LINE_SEARCH;
      }
      else
      {
        options_ &= ~LINE_SEARCH;
      }
    }

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

    bounds_->configure ( myProps );
  }

  if ( rundat_ )
  {
    rundat_->solver->configure ( props );
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void BFGSModule::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );


  myConf.set ( PropNames::MAX_ITER,      maxIter_    );
  myConf.set ( XProps::REFORM_ITER,      reformIter_ );
  myConf.set ( PropNames::TINY,          tiny_       );
  myConf.set ( PropNames::PRECISION,     precision_  );
  myConf.set ( XProps::LINE_SEARCH_TOL,  lsearchTol_ );
  myConf.set ( PropNames::MAX_INCR,      maxIncr_    );
  myConf.set ( XProps::MAX_RES_INCR,     maxResIncr_ );

  myConf.set ( PropNames::LINE_SEARCH,
               ((options_ & LINE_SEARCH) != 0)  );

  myConf.set ( PropNames::DELTA_CONS,
               ((options_ & DELTA_CONS)  != 0)  );

  bounds_->getConfig ( myConf );

  if ( rundat_ )
  {
    rundat_->solver->getConfig ( conf );
  }
}


//-----------------------------------------------------------------------
//   advance
//-----------------------------------------------------------------------


void BFGSModule::advance ( const Properties& globdat )
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


void BFGSModule::solve

  ( const Properties&  info,
    const Properties&  globdat )

{
  using jive::util::Globdat;

  if ( ! rundat_ )
  {
    notAliveError ( JEM_FUNC );
  }

  RunData_&  d = * rundat_;

  Work_      w ( *this, globdat );

  bool       result;


  print ( System::info( myName_ ),
          "Starting the non-linear solver `",
          myName_, "\' ...\n" );

  d.getExtVector    ( w.fext, globdat );
  w.updateBounds    ( globdat );
  w.applyBounds     ();
  w.updateMatrix    ( globdat );
  w.updateRscale    ( globdat );
  w.calcIncrement   ( maxIncr_ );
  w.zeroConstraints ( *d.cons  );

  result = solve_   ( w, globdat );

  w.restoreConstraints ( *d.cons );

  info.set ( SolverInfo::CONVERGED,  result  );
  info.set ( SolverInfo::ITER_COUNT, w.iiter );
  info.set ( SolverInfo::RESIDUAL,   w.rnorm );

  if ( Globdat::hasVariable( myName_, globdat ) )
  {
    Properties  vars = Globdat::getVariables ( myName_, globdat );

    vars.set ( SolverInfo::CONVERGED,  result  );
    vars.set ( SolverInfo::ITER_COUNT, w.iiter );
    vars.set ( SolverInfo::RESIDUAL,   w.rnorm );
  }

  if ( ! result )
  {
    throw Exception (
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
          "The non-linear solver converged in ",
          w.iiter, " iterations\n\n" );
}


//-----------------------------------------------------------------------
//   cancel
//-----------------------------------------------------------------------


void BFGSModule::cancel ( const Properties& globdat )
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


bool BFGSModule::commit ( const Properties& globdat )
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


void BFGSModule::setPrecision ( double eps )
{
  JEM_PRECHECK ( eps > 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double BFGSModule::getPrecision () const
{
  return precision_;
}


//-----------------------------------------------------------------------
//   setOption
//-----------------------------------------------------------------------


void BFGSModule::setOption

  ( Option  option,
    bool    yesno )

{
  options_.set ( option, yesno );
}


//-----------------------------------------------------------------------
//   setOptions
//-----------------------------------------------------------------------


void BFGSModule::setOptions ( Options options )
{
  options_ = options;
}


//-----------------------------------------------------------------------
//   getOptions
//-----------------------------------------------------------------------


BFGSModule::Options BFGSModule::getOptions () const
{
  return options_;
}


//-----------------------------------------------------------------------
//   setReformIter
//-----------------------------------------------------------------------


void BFGSModule::setReformIter ( idx_t reformIter )
{
  JEM_PRECHECK ( reformIter > -1 );

  reformIter_ = reformIter;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


idx_t BFGSModule::getReformIter () const
{
  return reformIter_;
}


//-----------------------------------------------------------------------
//   setMaxIter
//-----------------------------------------------------------------------


void BFGSModule::setMaxIter ( idx_t maxIter )
{
  JEM_PRECHECK ( maxIter > 0 );

  maxIter_ = maxIter;
}


//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


idx_t BFGSModule::getMaxIter () const
{
  return maxIter_;
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Module> BFGSModule::makeNew

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


void BFGSModule::declare ()
{
  using jive::app::ModuleFactory;

  ModuleFactory::declare ( TYPE_NAME,  & makeNew );
  ModuleFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   solve_
//-----------------------------------------------------------------------


bool BFGSModule::solve_

  ( Work_&             w,
    const Properties&  globdat )

{
  using std::fabs;
  using std::sqrt;

  RunData_&  d          = w.rundat;
  bool       lsearch    = false;
  int        smallSteps = 0;

  double     rmin;
  double     rmax;


  d.validMatrix = false;

  axpy           ( w.u, 1.0, w.du );
  w.applyBounds  ();
  d.updateModel  ( globdat );
  w.updateMatrix ( globdat );
  w.updateRscale ( globdat );

  print ( System::info( myName_ ), d.context,
          " : residual scale factor = ",
          d.nformat.print( w.rscale ), endl );

  if ( w.rscale <= tiny_ )
  {
    return true;
  }

  w.updateResidual ();

  rmin = rmax = w.rnorm;
  w.u0 = w.u;

  while ( true )
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

    if ( w.rnorm > 1.0e4 )
    {
      if ( options_ & LINE_SEARCH )
      {
        // Undo the last solution increment and execute the line
        // search procedure.

        d.validMatrix = false;
        w.u           = w.u0;

        w.applyBounds ();
        d.updateModel ( globdat );
        lineSearch_   ( w, globdat );
      }
      else
      {
        return false;
      }
    }

    w.calcIncrement ( maxIncr_ );

    if ( lsearch )
    {
      double  rnorm0 = w.rnorm;

      lineSearch_ ( w, globdat );

      /*if ( w.hasBounds() && (w.rnorm > (lsearchTol_ * rnorm0)) )
      {

        // Try again with a steepest descent search direction.

          lineSearch2_ ( w, globdat );

      }*/
    }
    else
    {
      w.uscale = 1.0;

      axpy             ( w.u, 1.0, w.du );
      w.applyBounds    ();
      d.updateModel    ( globdat );
      w.updateMatrix   ( globdat );
      w.updateResidual ();

      if ( (options_ & LINE_SEARCH) && (w.iiter >= 2) )
      {
        lsearch = true;
      }
    }

    if      (  d.validMatrix                   &&
              (w.jiter <  reformIter_)         &&
              (w.rnorm < (maxResIncr_ * rmin)) &&
              (w.rnorm <  rmax) )
    {
      double  dudr;
      double  dur0;

      w.dr = w.r0 - w.r;
      dudr = d.vspace->product ( w.du, w.dr );
      dur0 = d.vspace->product ( w.du, w.r0 );

      if ( (w.uscale * dudr * dur0) > 0.0 )
      {
        double  alpha = sqrt ( w.uscale * dudr / dur0 );

        w.pvecs[w.jiter] = w.dr + alpha * w.r0;
        w.qvecs[w.jiter] = (1.0 / dudr) * w.du;

        w.jiter++;
      }
    }
    else if ( reformIter_ > 0 )
    {
      d.validMatrix = false;

      w.updateMatrix   ( globdat );
      w.updateResidual ();
    }

    if ( fabs( w.uscale ) < 1.0e-2 )
    {
      smallSteps++;
    }

    if ( isTiny( w.uscale ) || (smallSteps > 5) )
    {
      lsearch    = false;
      smallSteps = 0;
    }

    rmax = max ( rmax, w.rnorm );
    rmin = min ( rmin, w.rnorm );
  }

  return false;
}


//-----------------------------------------------------------------------
//   lineSearch_
//-----------------------------------------------------------------------


void BFGSModule::lineSearch_

  ( Work_&             w,
    const Properties&  globdat )

{
  using std::fabs;

  const int     MAX_ITER = 10;
  const double  tol      = lsearchTol_;

  RunData_&     d        = w.rundat;
  double        enorm0   = d.vspace->product ( w.r0, w.du );
  int           iiter    = 0;

  double        r, r0, r1, r2;
  double        s, s1, s2;
  double        sbest;
  double        rbest;
  double        rmax;
  double        sgn0;
  double        sgn;
  double        sq;


  sbest = 0.0;
  rbest = 0.0;
  sgn0  = enorm0 > 0.0 ? 1.0 : -1.0;
  rmax  = tol  * fabs ( enorm0 );
  r0    = sgn0 * fabs ( enorm0 );
  r1    = 1.0;
  s     = 1.0;
  s1    = 0.0;
  s2    = 1.0;

  // Do a normal Newton step and check whether enough progress is
  // made. If not, then proceed with the line search algorithm.

  axpy             ( w.u, s, w.du );
  w.applyBounds    ();
  d.updateModel    ( globdat );
  w.updateMatrix   ( globdat );
  w.updateResidual ();

  if ( (fabs( w.enorm ) < rmax) || isTiny( enorm0 ) )
  {
    w.uscale = s;

    return;
  }

  print ( System::info( myName_ ), w.rundat.context,
          " : starting line search ...\n" );

  if ( (w.enorm - r0) * sgn0 > 0.0)
  {
    sbest =  0.0;
    rbest =  1.0;
    s     = -1.0;

    w.doTotalUpdate ( s, globdat );

    if ( fabs( w.enorm ) < rmax )
    {
      goto success;
    }

    r = w.enorm / r0;

    if (fabs(r) < fabs(rbest))
    {
      sbest = s;
      rbest = r;
    }
  }

  if ( (w.enorm - r0) * sgn0 < 0.0 )
  {
    // The residual norm decreases with the initial increment. Find
    // the interval within which the minimum residual norm is located
    // and then apply the secant method to find the location of the
    // minimum residual norm. First locate the upper bound.

    sbest = 1.0;
    r2    = w.enorm / r0;
    rbest = r2;
    iiter = 0;

    while ( (r2 > 0.0) && (iiter < MAX_ITER) )
    {
      s *= 2.0;

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r2 = w.enorm / r0;

      if ( fabs( r2 ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r2;
      }

      iiter++;
    }

    if ( iiter >= MAX_ITER )
    {
      // An upper bound has not been found. Try to find a better
      // solution with a smaller increment.

      if ( sbest > 1.0 )
      {
        goto best;
      }

      s = 1.0;

      for ( iiter = 0; iiter < MAX_ITER; iiter++ )
      {
        s /= 2.0;

        w.doTotalUpdate ( s, globdat );

        if ( fabs( w.enorm ) < rmax )
        {
          goto success;
        }

        r = w.enorm / r0;

        if ( fabs( r ) < fabs( rbest ) )
        {
          sbest = s;
          rbest = r;
        }
      }

      goto best;
    }

    s2    = s;
    iiter = 0;

    // Use the secant method to find the optimal increment.

    do
    {
      s = 0.5 * (s1 + s2);

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }

      sgn = r1 > 0.0 ? 1.0 : -1.0;
      sq  = r * r - r1 * r2;

      if ( sq > 0.0 )
      {
        sq = sqrt( sq );
        s  = s + (s - s1) * sgn * r / sq;
      }
      else
      {
        s = (r1 * s2 - r2 * s1) / (r1 - r2);
      }

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }

      if ( r > 0.0 )
      {
        r1 = r;
        s1 = s;
      }
      else
      {
        r2 = r;
        s2 = s;
      }

      iiter++;
    }
    while ( (iiter < MAX_ITER) && ! isTiny( r1 - r2 ) );
  }
  else
  {
    // The residual norm does not decrease with the initial increment.
    // Try to find a better solution with a smaller increment.

    s = 1.0;

    for ( iiter = 0; iiter < MAX_ITER; iiter++ )
    {
      s /= 2.0;

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }
    }

    // Search in the other direction.

    s = -1.0;

    for  ( iiter = 0; iiter < MAX_ITER; iiter++ )
    {
      s /= 2.0;

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }
    }
  }

 best:

  if ( (fabs( sbest ) < Float::EPSILON) && w.jiter )
  {
    // No improvement seems possible. Force the tangent matrix to
    // be reformed.

    sbest         = 0.0;
    d.validMatrix = false;
  }

  w.uscale = s = sbest;

  axpy             ( w.u, w.u0, s, w.du );
  w.applyBounds    ();
  d.updateModel    ( globdat );
  w.updateMatrix   ( globdat );
  w.updateResidual ();

  goto finally;

 success:

  w.updateMatrix   ( globdat );

 finally:

  w.uscale = s;

  print ( System::info( myName_ ), w.rundat.context,
          " : scale factor = ", d.nformat.print( w.uscale ),
          endl );
}


//-----------------------------------------------------------------------
//   lineSearch2_
//-----------------------------------------------------------------------


void BFGSModule::lineSearch2_

  ( Work_&             w,
    const Properties&  globdat )

{
  RunData_&  d      =  w.rundat;
  Vector     u0     =  w.u0;

  double     minExp = -20.0;
  double     maxExp =  10.0;

  double     unorm, rnorm;
  double     rmin,  rmax;
  double     s, s0, dt;

  idx_t      n;


  print ( System::info( myName_ ), w.rundat.context,
          " : starting steepest descent line search ...\n" );

  // Roll back to the previous solution.

  w.doTotalUpdate ( 0.0, globdat );

  // Use the steepest descent direction as the search direction.
  // Scale the search vector so that its magnitude is comparable
  // with the current solution.

  w.du  = w.r;
  unorm = d.vspace->norm2 ( u0 );
  rnorm = d.vspace->norm2 ( w.r );

  if ( (unorm > 0.0) && (rnorm > 0.0) )
  {
    w.du *= unorm / rnorm;
  }

  rmin = w.rnorm;
  rmax = max ( 0.5 * w.rnorm, precision_ );
  s0   = 0.0;
  n    = max (  7_idx, 2 * w.iiter );
  n    = min ( 31_idx, n );
  dt   = (maxExp - minExp) / (double) (n - 1);

  for ( idx_t i = 0; i < n; i++ )
  {
    s = std::exp ( minExp + (double) i * dt );

    w.doTotalUpdate ( s, globdat );

    if ( w.rnorm < rmax )
    {
      goto success;
    }

    if ( w.rnorm < rmin )
    {
      rmin = w.rnorm;
      s0   = s;
    }
  }

  // Apply the best update that has been found.

  s = s0;

  axpy             ( w.u, w.u0, s, w.du );
  w.applyBounds    ();
  d.updateModel    ( globdat );
  w.updateMatrix   ( globdat );
  w.updateResidual ();

  goto finally;

 success:

  w.updateMatrix ( globdat );

 finally:

  w.uscale = s;

  print ( System::info( myName_ ), w.rundat.context,
          " : scale factor = ", d.nformat.print( s ),
          endl );
}


//-----------------------------------------------------------------------
//   lineSearch3_
//-----------------------------------------------------------------------


void BFGSModule::lineSearch3_

  ( Work_&             w,
    const Properties&  globdat )

{
  using std::fabs;

  const int     MAX_ITER = 10;
  const double  tol      = lsearchTol_;

  RunData_&     d        = w.rundat;
  double        enorm0   = d.vspace->product ( w.r0, w.du );
  int           iiter    = 0;

  double        r, r0, r1, r2;
  double        s, s1, s2;
  double        sbest;
  double        rbest;
  double        rmax;
  double        sgn0;
  double        sgn;
  double        sq;


  sbest = 0.0;
  rbest = 0.0;
  sgn0  = enorm0 > 0.0 ? 1.0 : -1.0;
  rmax  = tol  * fabs ( enorm0 );
  r0    = sgn0 * fabs ( enorm0 );
  r1    = 1.0;
  s     = 1.0;
  s1    = -1.0;
  s2    = 1.0;

  // Do a normal Newton step and check whether enough progress is
  // made. If not, then proceed with the line search algorithm.

  axpy             ( w.u, s, w.du );
  w.applyBounds    ();
  d.updateModel    ( globdat );
  w.updateMatrix   ( globdat );
  w.updateResidual ();

  if ( (fabs( w.enorm ) < rmax) || isTiny( enorm0 ) )
  {
    w.uscale = s;

    return;
  }

  print ( System::info( myName_ ), w.rundat.context,
          " : starting negative line search ...\n" );

  if ( (w.enorm - r0) * sgn0 > 0.0)
  {
    sbest =  0.0;
    rbest =  1.0;
    s     = -1.0;

    w.doTotalUpdate ( s, globdat );

    if ( fabs( w.enorm ) < rmax )
    {
      goto success;
    }

    r = w.enorm / r0;

    if (fabs(r) < fabs(rbest))
    {
      sbest = s;
      rbest = r;
    }
  }

  if ( (w.enorm - r0) * sgn0 < 0.0 )
  {
    // The residual norm decreases with the initial increment. Find
    // the interval within which the minimum residual norm is located
    // and then apply the secant method to find the location of the
    // minimum residual norm. First locate the upper bound.

    sbest = 1.0;
    r2    = w.enorm / r0;
    rbest = r2;
    iiter = 0;

    while ( (r2 > 0.0) && (iiter < MAX_ITER) )
    {
      s *= 2.0;

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r2 = w.enorm / r0;

      if ( fabs( r2 ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r2;
      }

      iiter++;
    }

    if ( iiter >= MAX_ITER )
    {
      // An upper bound has not been found. Try to find a better
      // solution with a smaller increment.

      if ( sbest > 1.0 )
      {
        goto best;
      }

      s = 1.0;

      for ( iiter = 0; iiter < MAX_ITER; iiter++ )
      {
        s /= 2.0;

        w.doTotalUpdate ( s, globdat );

        if ( fabs( w.enorm ) < rmax )
        {
          goto success;
        }

        r = w.enorm / r0;

        if ( fabs( r ) < fabs( rbest ) )
        {
          sbest = s;
          rbest = r;
        }
      }

      goto best;
    }

    s2    = s;
    iiter = 0;

    // Use the secant method to find the optimal increment.

    do
    {
      s = 0.5 * (s1 + s2);

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }

      sgn = r1 > 0.0 ? 1.0 : -1.0;
      sq  = r * r - r1 * r2;

      if ( sq > 0.0 )
      {
        sq = sqrt( sq );
        s  = s + (s - s1) * sgn * r / sq;
      }
      else
      {
        s = (r1 * s2 - r2 * s1) / (r1 - r2);
      }

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }

      if ( r > 0.0 )
      {
        r1 = r;
        s1 = s;
      }
      else
      {
        r2 = r;
        s2 = s;
      }

      iiter++;
    }
    while ( (iiter < MAX_ITER) && ! isTiny( r1 - r2 ) );
  }
  else
  {
    // The residual norm does not decrease with the initial increment.
    // Try to find a better solution with a smaller increment.

    s = 1.0;

    for ( iiter = 0; iiter < MAX_ITER; iiter++ )
    {
      s /= 2.0;

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }
    }

    // Search in the other direction.

    s = -1.0;

    for  ( iiter = 0; iiter < MAX_ITER; iiter++ )
    {
      s /= 2.0;

      w.doTotalUpdate ( s, globdat );

      if ( fabs( w.enorm ) < rmax )
      {
        goto success;
      }

      r = w.enorm / r0;

      if ( fabs( r ) < fabs( rbest ) )
      {
        sbest = s;
        rbest = r;
      }
    }
  }

 best:

  if ( (fabs( sbest ) < Float::EPSILON) && w.jiter )
  {
    // No improvement seems possible. Force the tangent matrix to
    // be reformed.

    sbest         = 0.0;
    d.validMatrix = false;
  }

  w.uscale = s = sbest;

  axpy             ( w.u, w.u0, s, w.du );
  w.applyBounds    ();
  d.updateModel    ( globdat );
  w.updateMatrix   ( globdat );
  w.updateResidual ();

  goto finally;

 success:

  w.updateMatrix   ( globdat );

 finally:

  w.uscale = s;

  print ( System::info( myName_ ), w.rundat.context,
          " : scale factor = ", d.nformat.print( w.uscale ),
          endl );
}

JIVE_END_PACKAGE( implict )
