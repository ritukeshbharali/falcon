/* Include jem and jive headers */

#include <jive/app/ChainModule.h>
#include <jive/app/ControlModule.h>
#include <jive/app/OutputModule.h>
#include <jive/app/ReportModule.h>
#include <jive/app/Application.h>
#include <jive/app/InfoModule.h>
#include <jive/app/SampleModule.h>
#include <jive/app/UserconfModule.h>
#include <jive/geom/declare.h>
#include <jive/model/declare.h>
#include <jive/femodel/declare.h>
#include <jive/fem/declare.h>
#include <jive/fem/InputModule.h>
#include <jive/fem/InitModule.h>
#include <jive/fem/ShapeModule.h>
#include <jive/implict/LinsolveModule.h>
#include <jive/implict/NonlinModule.h>
#include <jive/implict/ArclenModule.h>
#include <jive/implict/Park3Module.h>
#include <jive/algebra/declare.h>
#include <jive/implict/declare.h>
#include <jive/app/declare.h>
#include <jive/gl/declare.h>
#include <jive/gl/DisplayModule.h>
#include <jive/gl/FemViewModule.h>
#include <jive/gl/GraphModule.h>

/* Include MPI relevant modules */

#if defined(WITH_MPI)
#include <jive/fem/MPInputModule.h>
#include <jive/fem/PartitionModule.h>
#endif

/* Include Falcon models */

#include "constraints/FalconConstraintModels.h"   // Constraints
#include "fem/basic/FalconBasicModels.h"          // Basic
#include "fem/solidmech/FalconSolidMechModels.h"  // Solid Mechanics
#include "fem/poromech/FalconPoroMechModels.h"    // Poro Mechanics
#include "fem/biomech/FalconBioMechModels.h"      // Bio Mechanics
#include "io/models/FalconIOModels.h"             // Input-Output

/* Include Falcon modules */

#include "io/modules/FalconIOModules.h"           // Input-Output
#include "steppers/FalconStepperModules.h"        // Steppers


using namespace jem;

using jive::app::Application;
using jive::app::Module;
using jive::app::ChainModule;
using jive::app::OutputModule;
using jive::app::InfoModule;
using jive::app::ControlModule;
using jive::app::ReportModule;
using jive::app::SampleModule;
using jive::app::UserconfModule;
using jive::fem::InputModule;
using jive::fem::InitModule;
using jive::fem::ShapeModule;
using jive::implict::ArclenModule;
using jive::implict::NonlinModule;
using jive::implict::LinsolveModule;
using jive::implict::Park3Module;
using jive::gl::FemViewModule;
using jive::gl::GraphModule;
using jive::gl::DisplayModule;

#if defined(WITH_MPI)
using jive::fem::MPInputModule;
using jive::fem::PartitionModule;
#endif

//-----------------------------------------------------------------------
//   mainModule
//-----------------------------------------------------------------------


Ref<Module> mainModule ()
{
  Ref<ChainModule>    chain = newInstance<ChainModule> ();
  Ref<ControlModule>  ctrl;

  /*
   * Declare internal shapes, models and matrix builders. These
   * functions essentially store pointers to construction functions
   * that are called when Jive needs to create a shape, model or
   * matrix builder of a particular type.
   */

  jive::geom::   declareIShapes   ();
  jive::geom::   declareShapes    ();
  jive::model::  declareModels    ();
  jive::fem::    declareMBuilders ();
  jive::implict::declareModels    ();
  jive::algebra::declareMBuilders ();
  jive::femodel::declareModels    ();

  /* Declare Falcon Models */

  declareFalconConstraintModels   ();
  declareFalconBasicModels        ();
  declareFalconSolidMechModels    ();
  declareFalconPoroMechModels     ();
  declareFalconBioMechModels      ();
  declareFalconIOModels           ();

  /* Declare Falcon Modules */

  declareFalconIOModules          ();
  declareFalconStepperModules     ();

  /* Declare all dynamically added modules in *.pro file,
   * not known at compile time.
   */

  jive::app     ::declareModules ();
  jive::implict ::declareModules ();
  jive::gl      ::declareModules ();

  // --------------------------------------------------------------
  //  add default modules
  // --------------------------------------------------------------

  #if defined(WITH_MPI) 
  chain->pushBack ( newInstance<MPInputModule>   ( ) );  
  chain->pushBack ( newInstance<PartitionModule> ( ) );  
  #endif
  
  chain->pushBack ( newInstance<InputModule>     ( ) );  
  chain->pushBack ( newInstance<ShapeModule>     ( ) );
  chain->pushBack ( newInstance<InitModule>      ( ) );
  chain->pushBack ( newInstance<InfoModule>      ( ) );

  // --------------------------------------------------------------
  // add user-defined modules
  // --------------------------------------------------------------

  chain->pushBack ( newInstance<UserconfModule> ( "extraModules" ) );

  // --------------------------------------------------------------
  // add control and report modules
  // --------------------------------------------------------------

  // The ControlModule controls when the program stops. 

   ctrl = newInstance<ControlModule> (
      "control",
        NIL,
        ControlModule::FG_MODE
        );

  ctrl ->runWhile ( "i > 0" );
  chain->pushBack ( ctrl );

  // Wrap chain module in a ReportModule that prints some
  // overall information about the current calculation.

  return newInstance<ReportModule> ( "report", chain );

  // Only for Mac OS

 // return newInstance<DisplayModule> ( newInstance<ReportModule> ( "report", chain ));
}


//-----------------------------------------------------------------------
//   main
//-----------------------------------------------------------------------


int main ( int argc, char** argv )
{
  #if defined(WITH_MPI) 
  return Application::pexec ( argc, argv, & mainModule );
  #else
  return Application::exec ( argc, argv, & mainModule );
  #endif
}
