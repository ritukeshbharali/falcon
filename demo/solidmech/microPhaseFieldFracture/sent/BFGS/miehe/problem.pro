/* Input file for Single Edge Notched specimen under Tension (SENT) */


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
  runWhile = "i < 161";
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
// 'multi', within which we define other models.

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = ["bulk","cons","lodi"];

    bulk = "MicroPhaseFracture"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "MiehePhase";
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
      penalty         = 500.;

      keepOffDiags    = false;
      arcLenMode      = false;

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
  modules = ["solver","graph","lodi","view","vtk","sample"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.025;
    startIncr      = 1.e-4;
    minIncr        = 1.e-12;
    reduceStep     = 54;

    solver = "BFGS"
    {
      precision  = 1.e-4;
      maxIter    = 5000;
      reformIter = 1;
      lineSearch  = true;
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
        matrixOrdering  = "metis";
        matrixScaling = true;
        parFactorize  = true;
        
        //type = "GMRES";
        //precon.type="ILUd";
        //precon.reorder = true;
        //precon.maxFill = 3.00000;
        //precon.dropTol = 1.00000e-08;
        //precon.diagShift = 0.00000;
        //precon.zeroThreshold = 1.00000e-08;
        //precon.minSize = 0;
        //precon.quality = 1.00000;
        //precision = 1.0e-08;

        //type = "SkylineLU";
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

  /*  
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.0e-6;
    
      maxIter   = 100;
  
      solver =
      {
        //type = "GMRES";
        type = "SkylineLU";
        lenient = true;
        useThreads = true;
      };
  };
  */

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


  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
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
                 //"model.model.lodi.load[0]" 
               ];
  };     

};
