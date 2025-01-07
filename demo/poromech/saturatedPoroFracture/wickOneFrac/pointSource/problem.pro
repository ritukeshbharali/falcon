/* 
  Case    : Hydraulic fracturing of specimen with single crack
  Ref     : Fig 2 in DOI: 10.1007/s10596-015-9532-5 
  FE model: SaturatedPorousFracture
  Material: MiehePhase (Spectral split based phase-field fracture)
  Loading : Neumann (Fluid injection in the crack)
  Implicit: FlexArclen solver (Nonlin and TSArclen)
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
  fgMode   = true;
  pause    = 0.;
  runWhile = "t < 1.25";
};


// 'input' refers to the Input Module. 'file' containing
// mesh, initial solution, and constraints are provided.

input =
{
  file = "$(CASE_NAME).data";
};


// 'init' refers to the Initialization Module. 'reorder'
// is set to false. The solution 'vectors', current solution
// 'state' is set to 0.0. One can complete ignore this as these
// vectors are initialized with 0.0.

init  =
  {
    reorder = false;
    vectors = [ "state = 0.0"];

  };


// Jive's concept of model and modules are used here. Models 
// perform the actions requested by the modules. We define a
// 'model' of type 'Matrix', the matrix model is of 'type'
// 'FEM'. Multiple sub-models may contribute to the system (
// stiffness) matrix and forces, so we choose 'model' as
// 'multi', within which we define models 'bulk', 'force1',
// and 'bc'. 'bulk' is of type 'SaturatedPorousFracture' and 
// operates on 'DomainElems'. 'shape' functions and integration
// scheme is then defined, along with other model parameters. 
// 'force' uses a LoadScale model defined in Jive in 
// conjunction with a Neumann model. 'bc' is ConstLoadArclen
// model used to scale the time-step size instead of an external
// load vector. In this model, the Neumann load is never scaled.

model = "Matrix"
{
  matrix.type = "FEM";
  // matrix.symmetric = true;

  model       =  "Multi"
  {
    models = [ "bulk", "force1", "bc"];

    bulk = "SaturatedPorousFracture"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };

      // Unit [N,m,s]
      
      material =
      {
        type   = "MiehePhase";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 1.e+09;    // 20,000 MPa
        poisson  = 0.2;
      };

      intrinPerm       = 1.e-12;
      fluidVisc        = 0.001;
      fracIntrinPerm   = 1.e-09;  // used as upper bounds
      fracFluidVisc    = 0.001;   // used as upper bounds
      solidStiff       = 1.e+09;
      fluidStiff       = 4.e+07;
      porosity         = 0.3;
      biotCoeff        = 1.0;
      dtime            = 0.001;
      gravity          = false;
      keepOffDiags     = true;
      arcLenMode       = true;


      permType         = "cubicIsoPerm";
      boundPerm        = true;
      randPorosity     = false;
      fracPorosity     = false;

      rhoSolid         = 2000.;
      rhoFluid         = 1000.;

      fractureType     = "at2";
      griffithEnergy   = 1.;
      lengthScale      = 0.05;
      tensileStrength  = 0.0;

    };

    force1 =  "LoadScale"
    {
       model =
       {
         type     = "Neumann";
         elements = "InjectionElems";
         loads    = [0.0,0.0,1.e-1,0.0];
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
    };

    bc =
    {
      type     = "ConstLoadArclen";
      optIter  = 6; //25; //100;
      maxIncr  = 0.1; // 5.0; // 1.e+1; // 16.e+0;
      minIncr  = 1.e-4; // 1.e-4; // 1.e+0;

      swtEnergy = 1.e-7;   // 1.e-4;
      swtIter   = 5;

    };

  };
};


// If you look into the chain module in main.cpp, the last 
// module is 'extraModules'. 'extraModules' allow the user to
// define a chain of modules to be executed at runtime. Here,
// we have added 'solver', 'view', 'vtk' and 'sample' modules.
// These modules would be executed in the order presented. 
// 'solver' module is of type 'FlexArclen', which adaptively 
// switches between Nonlin (Newton-Raphson solver) and 
// tsArcLen, which is the time-step size computing arc-length
// method. For both solvers, the inner linear 'solver', 
// is set to 'SkylineLU' from Jive. 'view' is of type
// FemView, which allows real-time visualization of FE
// solution fields and internal variables. 'vtk' is of type
// Paraview, which allows exporting solution and other fields
// into vtu file for visualization in Paraview, VisIt. 'sample'
// of type Sample module allows to export parameters such as the
// step number (i), time (t) and other variables to text files.

extraModules =
{
  modules = ["solver","view","vtk","sample"];
  
  solver  =
    {
      type  = "FlexArclen";
      writeStats    = true;
      timeStepBased = true;
      
      nonLin  =
      {
        maxIter = 25;
        precision = 1.e-4;
        reformIter = 0;
        lineSearch  = false;
        bounds  = ["b1"];
        solver  =
        {
          type        = "Pardiso";
          numThreads  = 2;
          sortColumns = 1;
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
        maxIter   = 25;
        precision = 1.e-4;
        loadScale = 0.001;
        solver  =
        {
          type        = "Pardiso";
          numThreads  = 2;
          sortColumns = 1;
        };
        model = none;
      };

      tsArcLen  =
      {
        deltaCons = false;
        maxIter   = 25;
        precision = 1.e-4;
        loadScale = 0.001;
        solver  =
        {
          type        = "Pardiso";
          numThreads  = 2;
          sortColumns = 1;
        };
        model = none;
      };
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
       fileName      = "$(CASE_NAME)_out";
       elements      = "DomainElems";
       printInterval = 1;
    };


  // Sample Module is used to export data stored in the global
  // database. Here we export step number (i), time(t), and the
  // nonlinear iterations required to converge (iterCount) for
  // every step upon acceptance of the solution.

  sample = "Sample"
  {
    // Save iteration data in file.
    file = "$(CASE_NAME)_iter.dat";
    dataSets = [ 
                 "i",
                 "t", 
                 "solverInfo.iterCount"
               ];
    sampleWhen = "accepted";           
  }; 
  

};
