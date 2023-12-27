/* 
  pm02    : Liakopoulos experiment
  FE model: Two Phase Unsaturated Porous
  Material: Hooke
  Loading : Neumann
  Implicit: Nonlin
  Solver  : Skyline
*/

log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 125";
};

input =
{
  file = "$(CASE_NAME).data";
};

init  =
  {
    reorder = false;
    vectors = [ "state = 0.0", "oldState = 0.0"];

  };

model = "Matrix"
{
  matrix.type = "FEM";
  // matrix.symmetric = true;

  model       =  "Multi"
  {
    models = [ "bulk", "force" ];

    bulk = "TwoPhaseUnsaturatedPorous"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Quad8";
        intScheme = "Gauss3*Gauss3";
      };

      // Unit [N,m,s]
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRESS";

        young    = 1.3e+6;
        poisson  = 0.4;
        rho      = 0.0;
      };

      retention.type = "Liakopoulos";

      intrinPerm    = 4.5e-13;
      fluidVisc     = 0.001;
      solidStiff    = 1.0e+10;
      fluidStiff    = 2.2e+09;
      porosity      = 0.2975;
      biotCoeff     = 1.0;
      dtime         = 60.0;

      rho_solid     = 2000.;
      rho_fluid     = 1000.;

    };

    force =  "LoadScale"
    {
       scaleFunc = 
       "
       save
       f_max   = 1.0
       ,
       t_max   = 5
       let
       t       = i
       return
       if ( t < t_max )
       f_max
       else
       0.0
       endif
       ";

       model =
       {
         type     = "Neumann";
         elements = "TopElems";
         loads    = [0.0,0.0,-1.0e-4];
         shape  =
          {
            type  = "BLine3";
            shapeFuncs  =
            {
              type  = "Quadratic";
            };
            intScheme = "Gauss2";
          };
       };
    };  

  };
};

extraModules =
{
  modules = ["solver","view"];
  
  solver = 
  {
      type      = "Nonlin";
      precision = 1.0e-4;    
      maxIter   = 5;

      solver =
      { 
        type        = "Skyline";
        useThreads  = true;
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
        data = "state[dp]";
      };
    };
    
    //updateWhen = "accepted";
    
  };

  vtk = "vtkWriter"
    {
       fileName   = "$(CASE_NAME)_out";
       elements = "DomainElems";
       interval = 1;
       data     = ["stress"];
       dataType = "nodes";

    };

};
