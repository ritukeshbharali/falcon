/* Input file for gradient polycrystal FE simulation */


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
  runWhile = "i < 1001";
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
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = ["bulk1","bulk2","cons","lodi"];

    bulk1 = "GradientCrystalPlasticity"

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

        young    = 205.e+03;
        poisson  = 0.29;
      };

      rotation = 10.0;

      dtime  = 0.1;

      slips = [ "slip0", "slip1" ];

      slip0 = 
      {
        n         = 9.0;
        tstar     = 1.0e+2;
        tauY      = 50.;
        selfH     = 100.;
        plane     = [0.9397,0.3420,0.];
        direction = [-0.3420,0.9397,0.];
      };

      slip1 =
      {
        n         = 9.0;
        tstar     = 1.0e+2;
        tauY      = 50.;
        selfH     = 100.;
        plane     = [0.7660,0.6428,0.];
        direction = [-0.6428,0.7660,0.];
      };
    };

    bulk2 = "GradientCrystalPlasticity"

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

        young    = 205.e+03;
        poisson  = 0.29;
      };

      rotation = 0.0;

      dtime  = 0.1;

      slips = [ "slip0", "slip1" ]; // rotated pi/6

      slip0 = 
      {
        n         = 9.0;
        tstar     = 1.0e+2;
        tauY      = 50.;
        selfH     = 100.;
        plane     = [0.9397,0.3420,0.];
        direction = [-0.3420,0.9397,0.];
      };

      slip1 =
      {
        n         = 9.0;
        tstar     = 1.0e+2;
        tauY      = 50.;
        selfH     = 100.;
        plane     = [0.7660,0.6428,0.];
        direction = [-0.6428,0.7660,0.];
      };
    };

    /*    
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "TopNodes";
      dofs       = "dx";
      factors    = [ 1.0];
      stepSize   = 1.e-5;
    };
    */

    cons =
       {
         type     = "Neumann";
         elements = "TopElems";
         loadInit  = 0.0;
         loadIncr  = 0.1;
         dof       = "dx";
         shape  =
          {
            type  = "BLine2";
            shapeFuncs  =
            {
              type  = "Linear";
            };
            intScheme = "Gauss1";
          };
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
  modules = ["solver","lodi","graph","view","view2","vtkout"];
  
  solver = 
  {
      type = "Nonlin";
      lineSearch = true;
      precision = 1.e-08;
      maxIter   = 100;
  
      solver =
      {
        type             = "Pardiso";
        sortColumns      = true;
        numThreads       = 4;
        lenient          = false;
        precision        = 1.e-10;
        matrixChecker    = false;
        matrixOrdering   = "metis";
        matrixScaling    = false;
        parFactorize     = true; 
      };
  };

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  ux | fx ";
      dataSets   = ["model.model.lodi.disp[0]","model.model.lodi.load[0]"];
    };

  graph =
  {
      type = "Graph";
      dataSets = "loadDisp";
      loadDisp =
      {
        key = "Load-displacement curve";
        xData = "model.model.lodi.disp[0]";
        yData = "model.model.lodi.load[0]";
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
    dataSets = [ "state", "tau" ];
    
    // 'state' is the solution vector found during simulations
    state =
    {
      type   = "Vector";
      vector = "state";
    };

    tau =
    {
      type   = "Table";
      vector = "nodes/tau";
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
        data = "tau[tau0]";
      };
    };
    
    //updateWhen = "accepted";
    
  };

  // FemViewModule, provides a visiualisation of the mesh during simulation. 
  view2 = "FemView"
  {

    window =
    {
      height = 300;
      width  = 600;
    };

    // Define the data sets to be visualized.
    dataSets = [ "state", "tau"];
    
    // 'state' is the solution vector found during simulations
    state =
    {
      type   = "Vector";
      vector = "state";
    };

    tau =
    {
      type   = "Table";
      vector = "nodes/tau";
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
        data = "tau[tau1]";
      };
    };
    
    //updateWhen = "accepted";
    
  };

  vtkout = "Paraview"
    {
       fileName      = "$(CASE_NAME)_out";
       elements      = "elems2DElems";
       printInterval = 1;
       cellData      = ["slip","tau","wip"];
    };
};
