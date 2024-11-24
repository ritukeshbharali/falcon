/* Input file for Winkler L-shaped panel */


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
  runWhile = "i < 801";
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
        dim    = 2;
        state  = "PLANE_STRESS";

        young    = 20.e+3;
        poisson  = 0.18;
        rho      = 1.0;
      };

      fractureType    = "czmCornelissen";
      griffithEnergy  = 0.113;
      lengthScale     = 10.;
      tensileStrength = 2.4;
      penalty         = 500.;

      keepOffDiags    = false;
      arcLenMode      = false;

    };
    
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "LoadNodes";
      dofs       = "dy";
      factors    = [ 1.0];
      stepSize   = 1.e-3;

    };

    // responsible for load-displacement curves

    lodi =
    {
      type   = "Lodi";
      group  = "LoadNodes";
    }; 

  };
};

extraModules =
{
  modules = ["solver","graph","lodi","vtk","sample"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.5;
    startIncr      = 1.e-3;
    minIncr        = 1.e-8;
    reduceStep     = 10000;

    solver = "Nonlin"
    {
      precision  = 1.e-4;
      maxIter    = 5000;
      reformIter = 0;
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
        xData = "-model.model.lodi.disp[1]";
        yData = "-0.1*model.model.lodi.load[1]";
      };
   };

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
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
