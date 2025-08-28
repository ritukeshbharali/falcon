/* Input file for local polycrystal FE simulation */


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
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 5";
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
// 'multi', within which we define models 'bulk1' to 'bulk5'
// for the different grains, 'cons' for the constraints and a
// 'lodi' for record the load-displacement on TopNodes.

model = "Matrix"
{
  matrix.type = "Sparse";
  matrix.symmetric = false;

  model       =  "Multi"
  {
    models = ["bulk1","bulk2","bulk3","bulk4","bulk5","cons","lodi"];

    bulk1 = "LocalCrystalPlasticity"

    {
      elements = "Domain1Elems";
      ipNodes  = "Domain1IPNodes";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+09;
        poisson  = 0.3;
      };

      dtime  = 1.0;

      slips = [ "slip0", "slip1" ];

      slip0 = 
      {
        n         = 3.0;
        tstar     = 1.0e-4;
        tauY      = 200.e+06;
        plane     = [1.,1.,0.];
        direction = [-1.,1.,0.];
      };

      slip1 =
      {
        n         = 7.0;
        tstar     = 1.0e-4;
        tauY      = 250.e+06;
        plane     = [0.,1.,0.];
        direction = [1.,0.,0.];
      };
    };

    bulk2 = "LocalCrystalPlasticity"

    {
      elements = "Domain2Elems";
      ipNodes  = "Domain2IPNodes";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+09;
        poisson  = 0.3;
      };

      dtime  = 1.0;

      slips = [ "slip0", "slip1" ]; // rotated pi/6

      slip0 = 
      {
        n         = 3.0;
        tstar     = 1.0e-4;
        tauY      = 200.e+06;
        plane     = [0.3660254,1.3660254,0.];
        direction = [-1.3660254,0.3660254,0.];
      };

      slip1 =
      {
        n         = 7.0;
        tstar     = 1.0e-4;
        tauY      = 250.e+06;
        plane     = [-0.5,0.8660254,0.];
        direction = [0.8660254,0.5,0.];
      };
    };

    bulk3 = "LocalCrystalPlasticity"

    {
      elements = "Domain3Elems";
      ipNodes  = "Domain3IPNodes";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+09;
        poisson  = 0.3;
      };

      dtime  = 1.0;

      slips = [ "slip0", "slip1" ];

      slip0 = 
      {
        n         = 3.0;
        tstar     = 1.0e-4;
        tauY      = 200.e+06;
        plane     = [0.12325683, 1.40883205,0.];
        direction = [-1.40883205, 0.12325683,0.];
      };

      slip1 =
      {
        n         = 7.0;
        tstar     = 1.0e-4;
        tauY      = 250.e+06;
        plane     = [-0.64278761, 0.76604444,0.];
        direction = [0.76604444, 0.64278761,0.];
      };
    };

    bulk4 = "LocalCrystalPlasticity"

    {
      elements = "Domain4Elems";
      ipNodes  = "Domain4IPNodes";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+09;
        poisson  = 0.3;
      };

      dtime  = 1.0;

      slips = [ "slip0", "slip1" ];

      slip0 = 
      {
        n         = 3.0;
        tstar     = 1.0e-4;
        tauY      = 200.e+06;
        plane     = [-0.3660254, 1.3660254,0.];
        direction = [-1.3660254, -0.3660254,0.];
      };

      slip1 =
      {
        n         = 7.0;
        tstar     = 1.0e-4;
        tauY      = 250.e+06;
        plane     = [-0.8660254, 0.5,0.];
        direction = [0.5, 0.8660254,0.];
      };
    };

    bulk5 = "LocalCrystalPlasticity"

    {
      elements = "Domain5Elems";
      ipNodes  = "Domain5IPNodes";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+09;
        poisson  = 0.3;
      };

      dtime  = 1.0;

      slips = [ "slip0", "slip1" ];

      slip0 = 
      {
        n         = 3.0;
        tstar     = 1.0e-4;
        tauY      = 200.e+06;
        plane     = [-0.70710678, 1.22474487,0.];
        direction = [-1.22474487, -0.70710678,0.];
      };

      slip1 =
      {
        n         = 7.0;
        tstar     = 1.0e-4;
        tauY      = 250.e+06;
        plane     = [-0.96592583, 0.25881905,0.];
        direction = [0.25881905, 0.96592583,0.];
      };
    };
    
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "TopNodes";
      dofs       = "dy";
      factors    = [ 1.0];
      stepSize   = 1.e-5;

    };

    // responsible for load-displacement curves

    lodi =
    {
      type   = "Lodi";
      group  = "TopNodes";
    }; 

  };
};


// If you look into the chain module in main.cpp, the last 
// module is 'extraModules'. 'extraModules' allow the user to
// define a chain of modules to be executed at runtime. Here,
// we have added 'solver', 'view', 'vtk' and 'sample' modules.
// These modules would be executed in the order presented. 
// 'solver' module is of type 'FlexArclen', which adaptively 
// switches between Nonlin (Newton-Raphson solver) and 
// tsArcLen, which is the time-step size computing arc-length
// method. For both solvers, the inner linear 'solver', 
// is set to 'SkylineLU' from Jive. 'view' is of type
// FemView, which allows real-time visualization of FE
// solution fields and internal variables. 'vtk' is of type
// Paraview, which allows exporting solution and other fields
// into vtu file for visualization in Paraview, VisIt.

extraModules =
{
  modules = ["solver","lodi","graph","view","vtkout"];
  
  solver = 
  {
      type = "Nonlin";
      //lineSearch = true;
    
      precision = 1.e-10;
    
      maxIter   = 50;
  
      solver =
      {
        type       = "SkylineLU";
        lenient    = true;
        useThreads = true;
      };
  };


  // Sample module to write the load-displacement data to
  // file.

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
    };


  // Graph for real-time visualization of 2D plot such as 
  // load-displacement (do not use this on cluster runs). 

  graph =
  {
      type = "Graph";
      dataSets = "loadDisp";
      loadDisp =
      {
        key = "Load-displacement curve";
        xData = "model.model.lodi.disp[1]";
        yData = "model.model.lodi.load[1]";
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
    dataSets = [ "state", "slip"];
    
    // 'state' is the solution vector found during simulations
    state =
    {
      type   = "Vector";
      vector = "state";
    };

    slip =
    {
      type   = "Table";
      vector = "nodes/slip";
    };
  
    // Settings for visualizing the finite element mesh.
    mesh =
    {
      // Apply a deformation to the mesh.
  
      deformation =
      {
        // Use the solution as the z-displacement.
        autoScale=false;
        scale = 0.01;
        dy = "state[dy]";
        dx = "state[dx]";
      };
  
      // Define the extra "plugins" for displaying more data.
  
      plugins = "colors";
  
      colors =
      {
        type = "MeshColorView";
        data = "slip[slip0]";
      };
    };
    
    //updateWhen = "accepted";
    
  };


  // Paraview Module is used to export solution and other
  // fields. By default, the solution field(s) are included.
  // In addition to that, one can export cellData (averaged
  // over Gauss point) and/or pointData (averaged over nodes).
  // One may also set a printInterval if one does not require
  // printing on every step. Here "allElems" is the set of
  // elements used for the export, which is usually the full
  // domain.

  vtkout = "Paraview"
    {
       fileName      = "$(CASE_NAME)_out";
       elements      = "allElems";
       printInterval = 1;
       cellData      = ["stress","strain","slip"];
       pointData     = ["stress","strain","slip"];
    };
};
