/* 
  cu01    : Winkler L-panel
  FE model: Linear Elasticity
  Material: Hooke
  Loading : Dirichlet
  Implicit: Nonlin
  Solver  : cuDSS
*/


// Setup log file for the entire simulation

log =
{
  pattern = "*.debug";
  file    = "$(CASE_NAME).log";
};


// Setup command line options and number of steps to run
// with runWhile, i being step number

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 5";
};


// Additional input file that stores constraints
// for the problem

input =
{
  file = "$(CASE_NAME).data";
};


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. In it, we define multi models.
   The "bulk" model is a Linear Elasticity FE model, which
   assembles stiffness matrix and internal force, among
   other things. The "cons" model of type Dirichlet enforces
   Dirichlet boundary conditions. The "lodi" model stores
   the load-displacement data for a certain set of nodes.
*/

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = ["bulk","cons","lodi"];

    bulk = "LinearElasticity"
    {
      elements = "DomainElems";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 20.e+3;
        poisson  = 0.18;
        rho      = 0.0;
      };
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
  modules = ["solver","lodi"];
  
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.e-6;
    
      maxIter   = 5;
  
      solver =
      {
        type    = "cuDSS";
      };
  };

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
    }; 
};
