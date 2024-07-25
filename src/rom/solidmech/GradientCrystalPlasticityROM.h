
/** @file GradientCrystalPlasticityROM.h
 *  @brief Crystal visco-plasticity reduced order model with gradient regularization.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 25 July 2024
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2024] (RB)
 * 
 *  @todo Add post-processing options (displacement, strains, stress)
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
#include <jem/util/Properties.h>
#include <jive/util/utilities.h>
#include <jive/util/XTable.h>
#include <jive/util/XDofSpace.h>
#include <jive/util/Printer.h>
#include <jive/util/Assignable.h>
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
using jem::util::SparseArray;
using jem::util::Properties;
using jive::util::XTable;
using jive::util::XDofSpace;
using jive::util::Assignable;
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

/** @brief Implements crystal visco-plasticity Reduced Order Model (ROM).
 * 
 *  The class \c GradientCrystalPlasticityROM implements a Reduced Order
 *  Model (ROM) for gradient-based crystal visco-plasticity simulations.
 *  The activity coefficients (corresponding to the chosen set of basis)
 *  are the unknown degrees of freedom, defined on a dummy modeNodes'
 *  NodeGroup. The corresponding finite element model class is 
 *  \c GradientCrystalPlasticityModel. 
 *  
 *  All parameters required to run the offline simulation with \c 
 *  GradientCrystalPlasticityModel are relevant for the online ROM
 *  \c GradientCrystalPlasticityROM. The user also needs to provide the
 *  basis/sensitivities file name and a vector containing the chosen
 *  basis.
 * 
 *  \note The number of nodes in the dummy 'modeNodes' NodeGroup must 
 *  match the size of the chosen 'basis' vector. Also, matrix.type
 *  must be set to "Sparse", instead of "FEM", symmetry set to 'false'.
 * 
 *  \note The code only works for single point Gauss integration.
 * 
 *  Below is an example how a model with three slip system is defined:
 * 
 *  \code
 *  matrix.type = "Sparse";
 *  matrix.symmetric = false; 
 * 
 *  models = ["bulk"];

    bulk = "GradientCrystalPlasticityROM"

    {
      elements   = "DomainElems";        // Element group
      modeNodes  = "modeNodes";          // (Dummy) Node group
      nmodes     = 7;                    // # of basis 
      nslips     = 2;                    // # of slips

      shape =
      {
        type      = "Triangle3";         // Element type
        intScheme = "Gauss1";            // Quadrature rule
      };

      dtime    = 1.0;                    // Time step-size
      tstar    = 1.e+02;                 // Relaxation time
      tauY     = 50.0;                   // Yield stress
      n        = 9.0;                    // Overstress function exponent

      slipHat0File = "pod_gamHat0.basis";  // slip sensitivity w.r.t unit load
      slipHatFile  = "pod_gamHat.basis";   // slip sensitivity w.r.t tau basis
      tauHatFile   = "pod_tau.basis";      // tau basis

      basis         = [0,1,2,3,4,5,6];     // chosen set of basis

 *  \endcode
 * 
 *  For the basis chosen in the input code, the the slip 
 *  sensitivity files are parsed as: pod_gamHat.basis0, 
 *  pod_gamHat.basis1, ...., pod_gamHat.basis6. The same for tau, 
 *  pod_tau.basis0, pod_tau.basis1, ...., pod_tau.basis6. These
 *  files must be present in the same directory as *.pro input
 *  file.
 *  
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

  static const char*        SLIP_HAT0_FILE;
  static const char*        BASIS;
  static const char*        SLIP_HAT_FILE;
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

  void                      getElemTau_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      getElemSlip_

    ( XTable&                 table,
      const Vector&           weights,
      const Vector&           state   );

  void                      getError_

    ( const Vector&           state,
      const Vector&           state0  );

  void                      checkCommit_

    ( const Properties&       params );
    

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

  String                     slipHat0File_;
  IdxVector                  basis_;
  String                     slipHatFile_;
  String                     tauHatFile_;

  Ref<FileReader>            iFile_;

  Matrix                     slipHat0_;  // ip x nslip
  Cubix                      slipHat_;   // ip x nslip x nbasis
  Cubix                      tauHat_;    // ip x nslip x nbasis

  // Pre-computations during the INIT stage

  Matrix                     M_;
  Vector                     fext_;

  int                        step_;

  // Error

  double                    err_;
  double                    errTotal_;

};