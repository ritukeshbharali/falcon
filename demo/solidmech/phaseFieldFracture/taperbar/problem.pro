/* Input file for Tapered Bar under Tension (TBT) */


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
  runWhile = "i < 51";
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
  modules = ["solver","graph","lodi","view","vtk"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.005;
    startIncr      = 1.e-2;
    minIncr        = 1.e-8;
    reduceStep     = 27;

    solver = "Nonlin"
    {
      precision = 1.e-3;
      maxIter = 100;
      lineSearch  = false;
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
        data = "state[phi]";
      };
    };
    
    //updateWhen = "accepted";
    
  };

  vtk = "Paraview"
    {
       fileName      = "$(CASE_NAME)_out";
       elements      = "DomainElems";
       printInterval = 10;
       pointData     = ["stress"];
    };

};
