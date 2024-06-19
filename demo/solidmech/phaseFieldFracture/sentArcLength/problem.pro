/* Input file for Single Edge Notched specimen under Tension (SENT) */


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
  runWhile = "i < 265";
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
   (TopNodes in this case).
*/

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

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

        young    = 210.e+3;
        poisson  = 0.3;
        rho      = 1.0;
      };

      fractureType    = "at2";
      griffithEnergy  = 2.7;
      lengthScale     = 0.015;
      tensileStrength = 0.0;
      arcLenMode      = true;
      keepOffDiags    = true;

    };

    cons =
    {
      type     = "DispArclen";
      optIter  = 10;
      dispIncr = 1.e-4;
      initDisp = 0.0;
      maxIncr  = 5.e-2;
      minIncr  = 1.e-4;

      swtEnergy = 5.e-3;
      swtIter   = 10;

      constraints = 
      {
        nodeGroups = [ "TopNodes" ];
        dofs = [ "dy" ];
        loaded = 0;
      };
    };
    
    /*
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "TopNodes";
      dofs       = "dy";
      factors    = [ 1.0];
      stepSize   = 1.e-5;

    };
    */

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
  modules = ["solver","graph","lodi","view","vtk","sample"];
  
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
          numThreads  = 4;
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
          numThreads  = 4;
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
        xData = "model.model.lodi.disp[1]";
        yData = "model.model.lodi.load[1]";
      };
      sampleWhen = "accepted";
   };


  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
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
    
    //updateWhen = "accepted";
    
  };

  vtk = "paraview"
    {
       fileName   = "$(CASE_NAME)_out";
       elements = "DomainElems";
       interval = 10;
       data     = ["pf"];
       dataType = "elems";

    };

  // sampleModule
  sample = "Sample"
  {
    // Save iteration data in file.
    file = "$(CASE_NAME)_iter.dat";
    dataSets = [ 
                 "i", 
                 "solverInfo.iterCount"
               ];
    sampleWhen  = "accepted";           
  };    

};
