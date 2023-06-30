/* 
  pm01    : Terzaghi consolidation
  FE model: Saturated Porous
  Material: Hooke
  Loading : Neumann
  Implicit: Nonlin
  Solver  : Skyline
*/

log =
{
  pattern = "*.debug";
  file    = "$(CASE_NAME).log";
};

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 2";
};

input =
{
  file = "$(CASE_NAME).data";
};

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = [ "bulk",  "force"];

    bulk = "SaturatedPorous"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Quad8";
        intScheme = "Gauss2*Gauss2";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 1.0e+08;
        poisson  = 0.0;
        rho      = 1.0;
      };

      intrin_perm   = 1.0e-14;
      fluid_visc    = 0.0089;
      solid_stiff   = 1.0e+10;
      fluid_stiff   = 2.0e+09;
      porosity      = 0.375;
      biot_coeff    = 1.0;
      dtime         = 9.1225;
    };
    

    force =  "LoadScale"
    {
       model =
       {
         type     = "Neumann";
         elements = "TopElems";
         loads    = [0.0,-1.0e+04,0.0];
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
      precision = 1.0e-6;    
      maxIter   = 100;

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
        scale = 0.001;
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
};
