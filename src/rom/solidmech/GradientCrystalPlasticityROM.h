
/** @file GradientCrystalPlasticityROM.h
 *  @brief Crystal visco-plasticity reduced order model with gradient regularization.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: XX YYYY 2024
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2024] (RB)
 * 
 *  @todo Simplify some data structures, reduce loops.
 * 
 *   Add props (n, tstar, tauY)
 *   Split gammaHat and tauHat files as: filename + incr
 *
 */

/* Include c++ headers */

#include <vector>
#include <cmath>

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/io/FileReader.h>
#include <jem/io/PrintWriter.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/sparse/utilities.h>
#include <jem/util/Flex.h>
#include <jem/util/Properties.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/util/Assignable.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/MatrixBuilder.h>
#include <jive/model/Model.h>
#include <jive/model/ModelFactory.h>
#include <jive/model/StateVector.h>
#include <jive/geom/InternalShape.h>
#include <jive/geom/IShapeFactory.h>
#include <jive/fem/ElementGroup.h>
#include <jive/util/Globdat.h>
#include <jive/fem/NodeGroup.h>

/* Include falcon headers */

#include "util/Arrays.h"
#include "util/ShapeUtils.h"
#include "materials/HookeMaterial.h"

using namespace jem;

using jem::io::FileReader;
using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::Flex;
using jem::util::SparseArray;
using jem::util::Properties;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::algebra::MatrixBuilder;
using jive::model::Model;
using jive::geom::IShape;
using jive::fem::NodeSet;
using jive::fem::ElementSet;
using jive::fem::ElementGroup;
using jive::fem::NodeGroup;
using jive::util::Globdat;
using std::vector;

//=======================================================================
//   typedefs
//=======================================================================

// Some handy aliases to avoid some typing

typedef MatmulChain<double,1>  MChain1;
typedef MatmulChain<double,2>  MChain2;
typedef MatmulChain<double,3>  MChain3;
typedef MatmulChain<double,4>  MChain4;

typedef ElementSet             ElemSet;
typedef ElementGroup           ElemGroup;

//=======================================================================
//   class GradientCrystalPlasticityROM
//=======================================================================

/** @brief Implements crystal visco-plasticity Reduced Order model.
 * 
 *  The class \c GradientCrystalPlasticityROM implements a crystal
 *  visco-plasticity finite element model with gradient-based 
 *  regularization of the slip(s). The displacements (2D, 3D) and the
 *  slip(s) are nodal degrees of freedom. The Schmid stress(es) are 
 *  defined on a dummy nodegroup 'ipNodes'. The model allows an 
 *  arbitrary number of slip planes, defined by the user at runtime.
 * 
 *  \note The number of nodes in the dummy 'ipNodes' nodegroup must match
 *  the total number of integration points in the model. Also, matrix.type
 *  must be set to "Sparse", instead of "FEM".
 * 
 *  Below is an example how a model with three slip system is defined:
 * 
 *  \code
 *  matrix.type = "Sparse";
 *  matrix.symmetric = false; 
 * 
 *  models = ["bulk","cons","lodi"];

    bulk = "GradientCrystalPlasticityROM"

    {
      elements = "DomainElems";        // Element group
      ipNodes  = "DomainIPNodes";      // (Dummy) Node group     

      shape =
      {
        type      = "Quad4";           // Element type
        intScheme = "Gauss2*Gauss2";   // Quadrature rule
      };
      
      material =
      {
        type   = "Hooke";              // Base material
        rank   = 2;                    // Rank (2D in this case!)
        state  = "PLANE_STRAIN";

        young    = 100.e+09;           // Young's modulus
        poisson  = 0.3;                // Poisson ratio
        rho      = 0.0;                // Material density
      };

      rotation = 45.0;                 // Slip system rotation
      dtime    = 1.0;                  // Time step-size

      slips = [ "slip0", "slip1", "slip2" ];

      slip0 = 
      {
        n         = 2.0;               // Exponent
        tstar     = 1.0;               // Relaxation time
        tauY      = 200.e+06;          // Yield stress
        selfH     = 100.;              // Self hardening
        plane     = [1.,1.,0.];        // Slip plane
        direction = [-1.,1.,0.];       // Slip direction
      };

      slip1 =
      {
        n         = 5.0;
        tstar     = 1.0;
        tauY      = 250.e+06;
        selfH     = 100.;
        plane     = [0.,1.,0.];
        direction = [1.,0.,0.];
      };

      slip2 =
      {
        n         = 10.0;
        tstar     = 0.5;
        tauY      = 350.e+06;
        selfH     = 100.;
        plane     = [1.,1.,0.];
        direction = [1.,-1.,0.];
      };
    };
 *  \endcode 
 */

class GradientCrystalPlasticityROM : public Model
{
  
 public:

  typedef GradientCrystalPlasticityROM  Self;
  typedef Model                         Super;

  static const char*        DISP_NAMES[3];
  static const char*        SHAPE_PROP;
  static const char*        DTIME_PROP;

  static const char*        NMODES_PROP;
  static const char*        NSLIPS_PROP;
  static const char*        MODE_NODES_PROP;       

  static const char*        GAMMA_HAT0_FILE;
  static const char*        BASIS;
  static const char*        GAMMA_HAT_FILE;
  static const char*        TAU_HAT_FILE;
  
                             GradientCrystalPlasticityROM

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )  const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );


 protected:

  virtual                  ~GradientCrystalPlasticityROM  ();

 private:

  void                      init_  ();

  void                      getMatrix_

    ( Ref<MatrixBuilder>      mbuilder,
      const Vector&           force,
      const Vector&           state,
      const Vector&           state0 );

  void                      getMatrix2_

    ( MatrixBuilder&          mbuilder ); 

  bool                      getTable_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getOutputData_

    ( Ref<XTable>             table,
      const Vector&           weights,
      const String&           contents,
      const Vector&           state );         

  void                      getStress_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      getElemStress_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      getStrain_

    ( XTable&                 table,
      const Vector&           weights );

  void                      getElemStrain_

    ( XTable&                 table,
      const Vector&           weights );

  void                      getTau_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      getElemTau_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      getElemSlip_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      checkCommit_

    ( const Properties&       params );  

  /* Initialize the mapping between integration points
   * and material points 
   */

  void                      initializeIPMPMap_ ();


 private:

  // FE mesh data

  Assignable<ElemGroup>      egroup_;
  Assignable<ElemSet>        elems_;
  Assignable<NodeSet>        nodes_;
  int                        rank_;
  Ref<IShape>                shape_;

  // Slips, time step-size, modes

  int                        nslips_;
  double                     dtime_;
  double                     tstar_;  
  double                     tauY_;
  double                     n_;

  // Dummy nodes for the modes

  int                        nmodes_;  
  String                     modeNGroup_;
  Assignable<NodeGroup>      modeNodes_;

  // Dof

  Ref<XDofSpace>             dofs_;
  IdxVector                  dofTypes_;
  IdxVector                  modeTypes_;

  // Sensitivities

  String                     gammaHat0File_;
  IdxVector                  basis_;
  String                     gammaHatFile_;
  String                     tauHatFile_;

  Ref<FileReader>            iFile_;

  Matrix                     gammaHat0_;  // ip x nslip
  Cubix                      gammaHat_;   // ip x nslip x nbasis
  Cubix                      tauHat_;     // ip x nslip x nbasis

  // Pre-computations during the INIT stage

  Matrix                     M_;
  Vector                     fext_;

  // Mapping between integration points and material points
  // (check if required!)

  SparseArray <int, 2>      ipMpMap_;

};