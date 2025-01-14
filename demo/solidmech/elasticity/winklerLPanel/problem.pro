/* Input file for Winkler L-panel elastic FE simulation */


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
  runWhile = "i < 11";
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

// 'extraModules' are modules that has the flexibility to be 
// defined at runtime. Here, we have added a 'solver' of type
// Nonlin, 'lodi' of Sample that writes a load-displacement
// curve, and 'vtk' of type Paraview to output solution, 
// stress, etc to vtu files. These modules would be executed 
// in the order presented.

extraModules =
{
  modules = ["solver","lodi","vtk"];
  
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.e-6;
    
      maxIter   = 5;
  
      solver =
      {
        type       = "SkylineLU";
        lenient    = true;
        useThreads = true;
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
       fileName      = "$(CASE_NAME)_out";
       elements      = "DomainElems";
       printInterval = 1;
       cellData      = ["stress"];
    };
};
