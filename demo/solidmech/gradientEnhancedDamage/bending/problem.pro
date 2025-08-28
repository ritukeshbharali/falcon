/* Input file for 3pt bending gradient damage FE simulation */

// Setup log file for the entire simulation

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
  fgMode  = true;
  pause = 0.;
  runWhile = "i < 175";
};


// 'input' refers to the Input Module. 'file' containing
// mesh, initial solution, and constraints are provided.

input =
{
  // Specify the name of the input data file.
  file = "problem.data";
};


// Jive's concept of model and modules are used here. Models 
// perform the actions requested by the modules. We define a
// 'model' of type 'Matrix', the matrix model is of 'type'
// 'FEM'. Multiple sub-models may contribute to the system (
// stiffness) matrix and forces, so we choose 'model' as
// 'multi', within which we define other models.

model =  "Matrix"
{
  matrix.type = "Sparse";
  model       = "Multi"
  {
    models = [ "all", "load", "lodi"];

    all = "GradientDamage"
    {
      elements  = "bulk";
      // thickness = 50.;

      // Material parameters.

      material = 
      {
	       type        = "Damage";
         rank        = 2;
	       state       = "PLANE_STRESS";   
         young       = 40e3;
         poisson     = 0.2;
         softening   = "exponential1";
	       eqvStrain   = "Mazars";
	       kappaI      = 0.000075;
         alpha       = 0.92;
         beta        = 300.;
         eta         = 10.; // von Mises prop
         lengthScale = 1.0;
      };
      
      // The shape object to be used. This depends on the finite
      // element mesh.

      shape =
      {
	      type      = "Triangle3";
        intScheme = "Gauss1";
      };

    };
 
    load =  "LoadScale"
    {
        scaleFunc = "1e-3 * (i)";
  
        model =
        {
          type     = "Constraints";
          conTable = "load";
        };
    };

    // responsible for load-displacement curves

    lodi =
    {
      type   = "Lodi";
      group  = "disp";
    };
  };
};


// 'extraModules' are modules that has the flexibility to be 
// defined at runtime. Here, we have added a 'solver' of type
// Nonlin, 'graph' of type Graph for real-time visualization 
// of load-displacement curve, 'lodi' of Sample that writes 
// the load-displacement curve, 'view' of type FemView for 
// real-time visualization of the damage and 'vtk' of type 
// Paraview to output solution, stress, etc to vtu files. 
// These modules would be executed in the order presented.

extraModules =
{
  modules = ["solver", "graph", "lodi", "view", "vtk" ];

  solver =  "Nonlin"
  {
     precision = 1.0e-5;
     maxIter   = 100;
     lineSearch = false;
     solver  = "SparseLU"
     {
         // useThreads = true;
     };
  };

  graph =
  {
      type = "Graph";
      dataSets = "loadDisp";
      loadDisp =
      {
        key = "Load-displacement curve";
        xData = "-model.model.lodi.disp[1]";
        yData = "-0.05*model.model.lodi.load[1]";
      };
   };
    
   // if you want to write the load-displacement data to file
    lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["-model.model.lodi.disp[1]","-0.05* model.model.lodi.load[1]"];
    };

view = "FemView"
{
  window =
  {
     height = 300;
     width  = 900;
  };

  // Define some data sets.

  dataSets = [ "disp" , "damage" ];

  disp =
  {
    type   = "Vector";
    vector = "state";
   };

   damage =
   {
      type  = "Table";
      table = "nodes/damage";
   };

   mesh =
   {
       deformation = "0.1 * disp";
       plugins = "colors";
       colors = {
          type = "MeshColorView";
          data = "damage";
          palette   = "custom";
          autoScale = false;
      };
    };
};

vtk = "Paraview"
    {
       fileName   = "$(CASE_NAME)";
       elements = "bulk";
       interval = 10;
       dofs     = ["dx","dy"];
       //data     = ["damage"];
       data     = ["stress","strain","damage"];

    };
};
