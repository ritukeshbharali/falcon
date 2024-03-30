/* Input file for Tapered Bar under Tension (TBT) */


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
  fgMode   = true;
  pause    = 0.;
  runWhile = "i < 1";
};


// Additional input file that stores constraints
// for the problem

input =
{
  file = "$(CASE_NAME).data";
};


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. In it, we define multi models.
   The "bulk" model is a PhaseFieldDamage FE model, which
   assembles stiffness matrix and internal force, among
   other things. The "cons" model of type Dirichlet enforces
   Dirichlet boundary conditions. The "lodi" model stores
   the load-displacement data for a certain set of nodes 
   (RightNodes in this case).
*/

model = "Matrix"
{
  matrix.type = "FEM";
  // matrix.symmetric = true;

  model       =  "Multi"
  {
    models = ["bulk","cons","lodi"];

    bulk = "GradientEnhancedDamage"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "ContinuumDamage";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 40.e3;
        poisson  = 0.2;

        softening  = "exponential1";
        equistrain = "strainEnergy";
        kappaI     = 0.000075;
        alpha      = 0.92;
        beta       = 300.;
      };

      lengthScale     = 0.75;

    };
    
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "RightNodes";
      dofs       = "dx";
      factors    = [ 1.0];
      stepSize   = 1.e-2;

    };

    // responsible for load-displacement curves

    lodi =
    {
      type   = "Lodi";
      group  = "RightNodes";
    }; 

  };
};

extraModules =
{
  modules = ["solver","graph","lodi","view"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.005;
    startIncr      = 1.e-5;
    minIncr        = 1.e-8;
    reduceStep     = 1000000000;

    solver = "Nonlin"
    {
      precision = 1.e-3;
      maxIter = 100;
      lineSearch  = false;
      solver =
      {
        type = "SkylineLU";
        useThreads=true;
      };
    };

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


  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[0]","model.model.lodi.load[0]"];
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
        data = "state[dx]";
      };
    };
    
    //updateWhen = "accepted";
    
  };

  vtk = "paraview"
    {
       fileName   = "$(CASE_NAME)_out";
       elements = "DomainElems";
       interval = 10;
       data     = ["stress"];
       dataType = "nodes";

    };

};
