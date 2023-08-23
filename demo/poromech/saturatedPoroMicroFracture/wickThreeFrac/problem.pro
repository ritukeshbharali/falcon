log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  fgMode   = true;
  pause    = 0.;
  runWhile = "i < 1000";
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
    models = [ "bulk", "force1" ];

    bulk = "SaturatedPorousMicroFracture"
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

      intrinPerm    = 1.e-12;
      fluidVisc     = 0.001;
      solidStiff    = 1.e+09;
      fluidStiff    = 4.e+07;
      porosity      = 0.3;
      biotCoeff     = 1.0;
      dtime         = 0.001;
      gravity       = false;
      keepOffDiags  = true;
      arcLenMode    = false;

      permType      = "cubicIsoPerm";
      pressurePsi   = false;
      randPorosity  = false;
      fracPorosity  = false;

      rhoSolid      = 2000.;
      rhoFluid      = 1000.;

      fractureType    = "at2";
      griffithEnergy  = 1.;
      lengthScale     = 0.05;
      tensileStrength = 0.0;

    };

    force1 =  "LoadScale"
    {
       model =
       {
         type     = "Neumann";
         elements = "InternalElems";
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
  };
};

extraModules =
{
  modules = ["solver","view","vtk"];
  
  solver = "ReduceStepping"
  {
    reduction      = 0.5;
    startIncr      = 1.e-3;
    minIncr        = 1.e-12;
    reduceStep     = 5000000;

  solver  = "Nonlin"
      {
        maxIter = 5000;
        precision = 1.e-3;
        reformIter = 0;
        lineSearch  = false;
        bounds  = ["b1"];
        solver  =
        {
          type        = "Pardiso";
          mtype       = 11;
          numThreads  = 5;
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
       printInterval = 1;
       cellData     = ["pf"];
    };

};
