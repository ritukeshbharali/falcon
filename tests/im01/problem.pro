/* 
  im01    : Tapered bar in tension
  FE model: Phase Fracture
  Material: Bourdin Phase
  Loading : Dirichlet, Arc-length
  Implicit: Flex Arc-length (Nonlin, Arc-length)
  Solver  : SkylineLU
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
  runWhile = "i < 115";
};


// Additional input file that stores constraints
// for the problem

input =
{
  file = "$(CASE_NAME).data";
};


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. In it, we define multi models.
   The "bulk" model is a PhaseField FE model, which
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

    bulk = "PhaseFracture"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "BourdinPhase";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+0;
        poisson  = 0.0;
        rho      = 1.0;
      };

      fractureType    = "at2";
      griffithEnergy  = 1.0;
      lengthScale     = 0.25;
      tensileStrength = 0.0;
      keepOffDiags    = true;
      arcLenMode      = true;

    };

    cons =
    {
      type     = "DispArclen";
      optIter  = 5;
      dispIncr = 1.e-2;
      initDisp = 0.0;
      maxIncr  = 5.e-3;
      minIncr  = 1.e-6;

      swtEnergy = 5.e-3;
      swtIter   = 10;

      constraints = 
      {
        nodeGroups = [ "RightNodes" ];
        dofs = [ "dx" ];
        loaded = 0;
      };
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
  
  solver  =
    {
      type  = "FlexArclen";
      
      nonLin  =
      {
        maxIter = 10;
        precision = 1.e-4;
        lineSearch  = false;
        bounds  = ["b1"];
        solver  =
        {
          type        = "Pardiso";
          mtype       = 11;
          numThreads  = 2;
          msglvl      = 0;
          sortColumns = 1;
          matrixChecker = 0;
        };
        b1 = 
         {
         dofType    = "phi";
         lowerBound = 0.;
         upperBound = 1.;
         };
      };
      arcLen  =
      {
        deltaCons = false;
        maxIter = 50;
        precision = 1.e-4;
        //loadScale = 0.00000;
        solver  =
        {
          type        = "Pardiso";
          mtype       = 11;
          numThreads  = 2;
          msglvl      = 0;
          sortColumns = 1;
          matrixChecker = 0;
        };
        model = none;
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
      sampleWhen  = "accepted";
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
