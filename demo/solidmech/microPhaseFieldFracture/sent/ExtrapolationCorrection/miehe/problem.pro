/* 
  Case    : Single Edge Notched specimen under Tension (SENT)
  Ref     : DOI: 10.1007/s00466-023-02380-1
  FE model: MicroPhaseFractureExtIt
  Material: MiehePhase (Spectral split based phase-field fracture)
  Loading : Dirichlet
  Implicit: Nonlin (Jive Newton-Raphson)
  Solver  : SkylineLU
*/


// Setup a log file for the entire simulation. For additional
// default options used in the analysis, always check the log.

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
// stiffness) matrix, forces and constraints, so we choose 
// 'model' as type 'multi'.

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


// If you look into the chain module in main.cpp, the last 
// module is 'extraModules'. 'extraModules' allow the user to
// define a chain of modules to be executed at runtime.

extraModules =
{
  modules = ["solver","graph","lodi","view","vtk","sample"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.025;
    startIncr      = 1.e-4;
    minIncr        = 1.e-12;
    reduceStep     = 54;

    solver = "Nonlin"
    {
      precision  = 1.e-4;
      maxIter    = 5000;
      reformIter = 0;
      lineSearch  = true;
      bounds = ["b1"];
      solver =
      {
        type = "SkylineLU";
        useThreads=true;
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


  // Paraview Module is used to export solution and other
  // fields. By default, the solution field(s) are included.
  // In addition to that, one can export cellData (averaged
  // over Gauss point) and/or pointData (averaged over nodes).
  // One may also set a printInterval if one does not require
  // printing on every step. Here "DomainElems" is the set of
  // elements used for the export, which is usually the full
  // domain.

  vtk = "Paraview"
    {
       fileName   = "$(CASE_NAME)_out";
       elements   = "DomainElems";
       interval   = 1;
       cellData   = ["pf"];
    };


  // Sample Module is used to export data stored in the global
  // database. Here we export step number (i), and the
  // nonlinear iterations required to converge (iterCount) for
  // every step upon acceptance of the solution.

  sample = "Sample"
  {
    // Save load displacement data in file.
    file = "$(CASE_NAME)_iter.dat";
    dataSets = [ 
                 "i", 
                 "solverInfo.iterCount"
               ];
  };     

};
