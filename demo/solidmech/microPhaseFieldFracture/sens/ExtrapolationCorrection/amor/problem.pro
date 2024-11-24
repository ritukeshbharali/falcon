/* Input file for Single Edge Notched specimen under Shear (SENS) */


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
  runWhile = "i < 501";
};


// Additional input file that stores constraints
// for the problem

input =
{
  file = "$(CASE_NAME).data";
};


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. matrix.symmetric is false 
   for MicroPhaseFractureExt (Extrapolation)
   "bulk" model. The "cons" model of type Dirichlet enforces
   Dirichlet boundary conditions. The "lodi" model stores
   the load-displacement data for a certain set of nodes 
   (TopNodes in this case).
*/

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = false;

  model       =  "Multi"
  {
    models = ["bulk","cons","lodi"];

    bulk = "MicroPhaseFractureExtIt"
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

        young    = 210.e+3;
        poisson  = 0.3;
        rho      = 1.0;
      };

      //keepOffDiags    = false;
      //arcLenMode      = false;

      fractureType    = "at2";
      griffithEnergy  = 2.7;
      lengthScale     = 0.015;
      tensileStrength = 0.0;
      penalty         = 150.;

    };
    
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "TopNodes";
      dofs       = "dx";
      factors    = [ 1.0];
      stepSize   = 1.e-4;

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
  modules = ["solver","graph","lodi","view"];
  
  solver = "ReduceStepping"
  {
    reduction      = .1;
    startIncr      = 1.e-4;
    minIncr        = 1.e-10;
    reduceStep     = 85;

    solver = "Nonlin"
    {
      precision   = 1.e-4;
      maxIter     = 5000;
      lineSearch  = false;
      bounds = ["b1"];
      solver =
      {
        type  = "Pardiso";
        lenient = false;
        mtype = 11;
        numThreads  = 4;
        msglvl  = 0;
        sortColumns = true;
        matrixChecker = false;
        matrixOrdering  = "amd";
        matrixScaling = true;
        parFactorize  = true;

        //type = "GMRES";
        //precon.type="ILUd";
        //precon.reorder = true;
        //precon.maxFill = 4.00000;
        //precon.dropTol = 1.00000e-08;
        //precon.diagShift = 0.00000;
        //precon.zeroThreshold = 1.00000e-08;
        //precon.minSize = 0;
        //precon.quality = 1.00000;
        //precision = 1.0e-08; 
        //type = "SparseLU";
        //useThreads=true;
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

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  ux | fx ";
      dataSets   = ["model.model.lodi.disp[0]","model.model.lodi.load[0]"];
    }; 


  // FemViewModule, provides a visualisation of the mesh during simulation. Remove this module when running on the cluster!
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
        scale = 0.1;
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

  vtk = "Paraview"
    {
       fileName   = "$(CASE_NAME)_out";
       elements   = "DomainElems";
       interval   = 1;
       cellData   = ["pf"];
    };

  sample = "Sample"
  {
    // Save load displacement data in file.
    file = "$(CASE_NAME)_iter.dat";
    //header = "  0.00000000e+00   0.00000000e+00   0.00000000e+00";
    dataSets = [ 
                 "i", 
                 "solverInfo.iterCount"
               ];
  };   

};
