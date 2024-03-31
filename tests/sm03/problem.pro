/* 
  sm03    : Beam bending
  FE model: Gradient Damage
  Material: Damage
  Loading : Dirichlet
  Implicit: Nonlin
  Solver  : Jive Skyline
*/

log =
{
  // Print informational messages to the terminal.
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};


input =
{
  // Specify the name of the input data file.
  file = "problem.data";
};

control =
{
  fgMode  = false;
  pause = 0.;
  runWhile = "i < 100";
};


model =  "Matrix"
{
  matrix.type = "Sparse";
  model       = "Multi"
  {
    models = [ "all", "load" ];

    all = "GradientDamage"
    {
      elements  = "bulk";
      // thickness = 50.;

      // Material parameters.

      material = 
      {
	       type        = "Damage";
         rank        = 2;
	       state       = "PLANE_STRESS";   
         young       = 40.e3;
         poisson     = 0.2;
         softening   = "exponential1";
	       eqvStrain   = "Mazars";
	       kappaI      = 0.000075;
         alpha       = 0.92;
         beta        = 300.;
         eta         = 10.; // von Mises prop
         lengthScale = 1.0;
      };
      
      // The shape object to be used. This depends on the finite
      // element mesh.

      shape =
      {
	      type      = "Triangle3";
        intScheme = "Gauss1";
      };

    };
 
    load =  "LoadScale"
    {
        scaleFunc = "1.5e-3 * (i)";
  
        model =
        {
          type     = "Constraints";
          conTable = "load";
        };
    };

  };
};
 
extraModules =
{
  modules = ["solver", "view" ];

  solver =  "Nonlin"
  {
     precision = 1.0e-5;
     maxIter   = 100;
     lineSearch = false;
     solver  = "SkylineLU"
     {
       useThreads = true;
     };
  };

  view = "FemView"
  {
    window =
    {
      height = 300;
      width  = 900;
    };

    // Define some data sets.

    dataSets = [ "disp" , "damage" ];

    disp =
    {
      type   = "Vector";
      vector = "state";
     };

     damage =
     {
       type  = "Table";
       table = "nodes/damage";
     };

     mesh =
     {
        deformation = "0.1 * disp";
        plugins = "colors";
        colors = {
          type = "MeshColorView";
          data = "damage";
          palette   = "custom";
          autoScale = false;
        };
      };
  };
};
