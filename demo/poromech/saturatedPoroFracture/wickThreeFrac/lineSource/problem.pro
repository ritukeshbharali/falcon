/* 
  Case    : Hydraulic fracturing of specimen with single crack
  Ref     : Fig 2 in DOI: 10.1007/s10596-015-9532-5 
  FE model: SaturatedPorousFracture
  Material: MiehePhase (Spectral split based phase-field fracture)
  Loading : Neumann (Fluid injection in the crack)
  Implicit: FlexArclen solver (Nonlin and TSArclen)
  Solver  : SkylineLU
*/

log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  fgMode   = true;
  pause    = 0.;
  runWhile = "t < 1.25";
};

input =
{
  file = "$(CASE_NAME).data";
};

init  =
  {
    reorder = false;
    vectors = [ "state = 0.0"];

  };

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
         loads    = [0.0,0.0,1.e-2,0.0];
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
       printInterval = 1;
    };

  // sampleModule
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
