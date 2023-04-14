/* 
  sm03    : 2D square 
  FE model: Linear Elasticity
  Material: Hooke
  Loading : Dirichlet
  Implicit: Nonlin
  Solver  : Intel Pardiso
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
    models = [ "bulk",  "force", "lodi"];

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
        state  = "PLANE_STRESS";

        young    = 1.e+4;
        poisson  = 0.3;
      };

    };

    force = 
    {
      type       = "Dirichlet";
      nodeGroups = "TopNodes";
      dofs       = "dy";
      factors    = [ 1.0];
      stepSize   = 0.0001;

    };

    lodi =
    {
      type   = "Lodi";
      group  = "TopNodes";
    };    

  };
};

extraModules =
{
  modules = ["solver","graph","view","vtk"];
   
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.0e-6;
    
      maxIter   = 100;
  
      solver =
      {
        type  = "Pardiso";
        mtype = 2;
        numThreads = 1;
        msglvl = 1;
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
        scale = 0.1;
        dy = "state[dy]";
        dx = "state[dx]";
      };
  
      // Define the extra "plugins" for displaying more data.
  
      plugins = "colors";
  
      colors =
      {
        type = "MeshColorView";
        data = "state[dy]";
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
