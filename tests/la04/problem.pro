/* 
  la04    : Tapered bar in tension (Umfpack)
  FE model: Phase Fracture Ext
  Material: Amor Phase
  Loading : Dirichlet
  Implicit: Nonlin
  Solver  : Umfpack
*/

log =
{
  pattern = "*.debug";
  file    = "$(CASE_NAME).log";
};

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 51";
};

input =
{
  file = "$(CASE_NAME).data";
};

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = false;

  model       =  "Multi"
  {
    models = ["bulk","cons","lodi"];

    bulk = "PhaseFractureExt"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "AmorPhase";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+0;
        poisson  = 0.0;
      };

      fractureType    = "at2";
      griffithEnergy  = 1.0;
      lengthScale     = 0.25;
      tensileStrength = 0.0;
      penalty         = 6500.;

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
  modules = ["solver","graph","view"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.005;
    startIncr      = 1.e-2;
    minIncr        = 1.e-8;
    reduceStep     = 27;

    solver = "Nonlin"
    {
      precision = 1.e-4;
      maxIter = 100;
      lineSearch  = false;
      bounds = ["b1"];
      solver =
      {
        type = "Umfpack";
      };

      b1 = 
      {
         dofType    = "phi";
         lowerBound = 0.;
         upperBound = 1.;
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
        data = "state[phi]";
      };
    };
    
    updateWhen = "accepted";
    
  };
};
