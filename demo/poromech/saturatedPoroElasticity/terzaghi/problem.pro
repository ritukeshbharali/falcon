/* 
  pm01    : Terzaghi consolidation
  FE model: Saturated Porous
  Material: Hooke
  Loading : Neumann
  Implicit: Nonlin
  Solver  : Skyline
*/

// Setup a log file for the entire simulation. For additional
// default options used in the analysis, always check the log.

log =
{
  pattern = "*.debug";
  file    = "$(CASE_NAME).log";
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
  runWhile = "i < 2";
};


// 'input' refers to the Input Module. 'file' containing
// mesh, initial solution, and constraints are provided.

input =
{
  file = "$(CASE_NAME).data";
};


// Jive's concept of model and modules are used here. Models 
// perform the actions requested by the modules. We define a
// 'model' of type 'Matrix', the matrix model is of 'type'
// 'FEM'. Multiple sub-models may contribute to the system (
// stiffness) matrix and forces, so we choose 'model' as
// 'multi', within which we define model 'bulk' and 'force'.
// 'bulk' is of type 'SaturatedPorous' and operates on
// 'DomainElems'. 'shape' functions and integration scheme 
// is then defined, along with other model parameters. 'force'
// uses a LoadScale model defined in Jive in conjunction with
// a Neumann model.

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = [ "bulk",  "force"];

    bulk = "SaturatedPorous"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Quad8";
        intScheme = "Gauss2*Gauss2";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 1.0e+08;
        poisson  = 0.0;
        rho      = 1.0;
      };

      intrinPerm    = 1.0e-14;
      fluidVisc     = 0.0089;
      solidStiff    = 1.0e+10;
      fluidStiff    = 2.0e+09;
      porosity      = 0.375;
      biotCoeff     = 1.0;
      dtime         = 9.1225;
    };
    

    force =  "LoadScale"
    {
       model =
       {
         type     = "Neumann";
         elements = "TopElems";
         loadInit  = -1.0e+4;
         loadIncr  = 0.0;
         dof       = "dy";
         shape  =
          {
            type  = "BLine3";
            shapeFuncs  =
            {
              type  = "Quadratic";
            };
            intScheme = "Gauss2";
          };
       };
    };  

  };
};

// If you look into the chain module in main.cpp, the last 
// module is 'extraModules'. 'extraModules' allow the user to
// define a chain of modules to be executed at runtime. Here,
// we have added 'solver' and 'view' modules. These
// modules would be executed in the order presented. 'solver'
// modules is of type 'Nonlin', which is the Newton-Raphson
// solver from Jive. It requires an inner linear 'solver', 
// which is set to 'SkylineLU' from Jive. 'view' is of type
// FemView, which allows real-time visualization of FE
// solution fields and internal variables.

extraModules =
{
  modules = ["solver","view"];
  
  solver = 
  {
      type      = "Nonlin";
      precision = 1.0e-6;    
      maxIter   = 100;

      solver =
      { 
        type        = "Skyline";
        useThreads  = true;
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
      // Apply a deformation to the mesh.
  
      deformation =
      {
        // Use the solution as the z-displacement.
        autoScale=false;
        scale = 0.001;
        dy = "state[dy]";
        dx = "state[dx]";
      };
  
      // Define the extra "plugins" for displaying more data.
  
      plugins = "colors";
  
      colors =
      {
        type = "MeshColorView";
        data = "state[dp]";
      };
    };
    
    //updateWhen = "accepted";
    
  };
};
