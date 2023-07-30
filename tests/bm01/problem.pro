/* 
  bm01    : Prostate Cancer
  FE model: Prostate Cancer
  Material: None
  Loading : Body 
  Implicit: Nonlin
  Solver  : Jive Skyline
*/


// Setup log file for the entire simulation

log =
{
  pattern = "*.debug";
  file    = "$(CASE_NAME).log";
};


// Setup command line options and number of steps to run
// with runWhile, i being step number

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 5";
};


// Additional input file that stores constraints
// for the problem

input =
{
  file = "$(CASE_NAME).data";
};

init  =
  {
    reorder = false;
    vectors = [ "state = 0.0", "oldState = init"];

  };


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. In it, we define multi models.
   The "bulk" model is a ProstateCancer FE model, which
   assembles stiffness matrix and internal force, among
   other things.
*/

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

extraModules =
{
  modules = ["solver","view"];
    
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

  // FemViewModule, provides a visiualisation of the mesh during simulation. 
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

};
