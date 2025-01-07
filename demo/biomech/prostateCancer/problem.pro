/* Input file for Prostate Cancer FE Model demo */

// Setup a log file for the entire simulation. For additional
// default options used in the analysis, always check the log.

log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};


// 'control' refers to the Control Module. fgMode = false 
// indicates that once the runWhile is fulfilled, the
// simulation will terminate. If set to true, the user can
// execute 'step 100' to run 100 more steps. In 'runWhile',
// 'i' is the step number. One can also use 't' for time, if
// simulation involves time.

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 1000";
};

// 'input' refers to the Input Module. 'file' containing
// mesh, initial solution, and constraints are provided.

input =
{
  file = "$(CASE_NAME).data";
};

// 'init' refers to the Initialization Module. 'reorder'
// is set to false. The solution 'vectors', current solution
// 'state' is set to 0.0 and previous solution 'oldState' is
// set to 'init'. 'init' is a NodeTable located in the file
// 'problem.init'. Note that 'problem.init' is sourced through
// Input Module where file = problem.data is loaded.

init  =
  {
    reorder = false;
    vectors = [ "state = 0.0", "oldState = init"];

  };

// Jive's concept of model and modules are used here. Models 
// perform the actions requested by the modules. We define a
// 'model' of type 'Matrix', the matrix model is of 'type'
// 'FEM' and is 'symmetric'. Multiple sub-models may contribute
// to the system (stiffness) matrix, so we choose 'model' as
// 'multi', within which we define a model 'bulk'. 'bulk' is of
// type 'ProstateCancer' and operates on 'DomainElems'. 'shape'
// functions and integration scheme is then defined, along with
// other model parameters.

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = ["bulk"];

    bulk = "ProstateCancer"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };

      lambda              = 160000.;
      tau                 = 0.0100000;
      growthRate          = 600.000;
      apoptosisRate       = 600.000;
      epsilon             = 5.00000e+06;
      nutrientAvg         = 1003.75;
      nutrientDev         = 73.0000;
      nutrientConsumption = 1003.75;
      nutrientDecay       = 1000.00;
      dtime               = 0.001;
      generateNewSeed     = false;
      
    };
  };
};

// If you look into the chain module in main.cpp, the last 
// module is 'extraModules'. 'extraModules' allow the user to
// define a chain of modules to be executed at runtime. Here,
// we have added 'solver', 'view' and 'vtk' modules. These
// modules would be executed in the order presented. 'solver'
// modules is of type 'Nonlin', which is the Newton-Raphson
// solver from Jive. It requires an inner linear 'solver', 
// which is set to 'SkylineLU' from Jive. 'view' is of type
// FemView, which allows real-time visualization of FE
// solution fields and internal variables. 'vtk' is of type
// Paraview, which allows exporting solution and other fields
// into vtu file for visualization in Paraview, VisIt.

extraModules =
{
  modules = ["solver","view","vtk"];
    
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.e-4;
      maxIter   = 10;
  
      solver =
      {
        type = "SkylineLU";
        lenient = true;
        useThreads = true;
      };
  };

  // FemView for real-time visualization (do not use this on
  // cluster runs). 

  view = "FemView"
  {

    window =
    {
      height = 300;
      width  = 600;
    };

    // Define the data sets to be visualized.
    dataSets = [ "state"];
    
    // 'state' is the solution vector found during simulations
    state =
    {
      type   = "Vector";
      vector = "state";
    };
  
    // Settings for visualizing the finite element mesh.
    mesh =
    {  
      // Define the extra "plugins" for displaying more data.
  
      plugins = "colors";
  
      colors =
      {
        type = "MeshColorView";
        data = "state[phi]";
      };
    };

    //updateWhen = "accepted";
    
  };

  // Paraview Module is used to export solution and other
  // fields. By default, the solution field(s) are included.
  // In addition to that, one can export cellData (averaged
  // over Gauss point) and/or pointData (averaged over nodes).
  // One may also set a printInterval if one does not require
  // printing on every step. Here "DomainElems" is the set of
  // elements used for the export, which is usually the full
  // domain.

  vtk = "Paraview"
    {
       fileName      = "$(CASE_NAME)_out";
       elements      = "DomainElems";
       printInterval = 10;
       cellData      = ["source"];
    };  

};
