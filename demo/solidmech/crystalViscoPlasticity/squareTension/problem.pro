/* Input file for block */


// Setup log file for the entire simulation

log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};


// Setup command line options and number of steps to run
// with runWhile, i being step number

control = 
{
  fgMode   = false;
  pause    = 0.;
  runWhile = "i < 10";
};


// Additional input file that stores constraints
// for the problem

input =
{
  file = "$(CASE_NAME).data";
};


/* Model tree for the the problem. We work with 'Matrix'
   type model of type 'FEM'. In it, we define multi models.
   The "bulk" model is a PhaseFieldDamage FE model, which
   assembles stiffness matrix and internal force, among
   other things. The "cons" model of type Dirichlet enforces
   Dirichlet boundary conditions. The "lodi" model stores
   the load-displacement data for a certain set of nodes 
   (TopNodes in this case).
*/

model = "Matrix"
{
  matrix.type = "FEM";
  matrix.symmetric = true;

  model       =  "Multi"
  {
    models = ["bulk","cons","lodi"];

    bulk = "CrystalViscoPlasticity"

    {
      elements = "DomainElems";

      shape =
      {
        type      = "Quad4";
        intScheme = "Gauss2*Gauss2";
      };
      
      material =
      {
        type   = "Hooke";
        rank   = 2;
        state  = "PLANE_STRAIN";

        young    = 100.e+09;
        poisson  = 0.3;
        rho      = 0.0;
      };

      dtime  = 1.0;

      slips = [ "slip0", "slip1", "slip2" ];

      slip0 = 
      {
        n         = 2.0;
        tstar     = 1.0;
        tauY      = 200.e+06;
        plane     = [1.,1.,0.];
        direction = [-1.,1.,0.];
      };

      slip1 =
      {
        n         = 5.0;
        tstar     = 1.0;
        tauY      = 250.e+06;
        plane     = [0.,1.,0.];
        direction = [1.,0.,0.];
      };

      slip2 =
      {
        n         = 10.0;
        tstar     = 0.5;
        tauY      = 350.e+06;
        plane     = [1.,1.,0.];
        direction = [1.,-1.,0.];
      };
    };
    
    cons = 
    {
      type       = "Dirichlet";
      nodeGroups = "TopNodes";
      dofs       = "dy";
      factors    = [ 1.0];
      stepSize   = 1.e-4;

    };

    // responsible for load-displacement curves

    lodi =
    {
      type   = "Lodi";
      group  = "TopNodes";
    }; 

  };
};

extraModules =
{
  modules = ["solver","lodi","view", "vtk"];
  
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.e-3;
    
      maxIter   = 50;
  
      solver =
      {
        type       = "Pardiso";
        //numThreads = 4;
        //lenient    = true;
        useThreads = true;
        //sortColumns = 1;
      };

      /*
      solver = "GMRES"
      {
        precon = "Dual"
        {
          precon1 = "ILUd"
          {
            reorder = true;
          };
          precon0 = "Coarse"
          {
            restrictor = "RigidBody"
            {
            };
          };
        };
        printInterval = 0.100000;
        noiseLevel    = 2;
        lenient       = true;
        precision     = 1.e-08;
        maxIter       = 2000;
        updatePolicy  = "Auto";
        restartIter   = 200;
      };
      */
  };

  lodi = "Sample"
    {
      file       = "$(CASE_NAME)_lodi.dat";
      header     = "  uy | fy ";
      dataSets   = ["model.model.lodi.disp[1]","model.model.lodi.load[1]"];
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
    dataSets = [ "state", "tau"];
    
    // 'state' is the solution vector found during simulations
    state =
    {
      type   = "Vector";
      vector = "state";
    };

    tau =
    {
      type   = "Table";
      vector = "nodes/tau";
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
        data = "tau[tau0]";
      };
    };
    
    //updateWhen = "accepted";
    
  };

  vtk = "paraview"
    {
       fileName      = "$(CASE_NAME)_out";
       elements      = "DomainElems";
       printInterval = 1;
    };
};
