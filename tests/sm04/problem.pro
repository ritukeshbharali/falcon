/* 
  sm05    : Block under tension
  FE model: Local Crystal Plasticity
  Material: Hooke (elastic part)
  Loading : Dirichlet
  Implicit: Nonlin
  Solver  : Jive Skyline
*/


// Setup log file for the entire simulation

log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
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


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. In it, we define multi models.
   The "bulk" model is a GradientCrystalPlasticity FE model, 
   which assembles stiffness matrix and internal force, among
   other things. The "cons" model of type Dirichlet enforces
   Dirichlet boundary conditions. The "lodi" model stores
   the load-displacement data for a certain set of nodes 
   (TopNodes in this case).
*/

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

extraModules =
{
  modules = ["solver","lodi","graph","view"];
  
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

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
    };

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

  // FemViewModule, provides a visiualisation of the mesh during simulation. 
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
};
