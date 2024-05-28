/** @file AbaqusMaterial.h
 *  @brief Abaqus UMAT interface.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 28 May 2024
 *  
 *  Updates (when, what and who)
 *     - [XX YYYY 2024]
 * 
 *  @warning Not thoroughly tested. Not all UMAT
 * function arguments are used.
 */

#ifndef ABAQUSMATERIAL_H
#define ABAQUSMATERIAL_H

#include <jive/Array.h>
#include <jem/base/String.h>
#include <jem/util/Flex.h>
#include <vector>

#include "Material.h"

using jem::String;
using jem::util::Flex;
using jive::Vector;
using jive::IntVector;
using std::vector;

// Define the Abaqus UMAT function signature
typedef void (*umatFunc)(
    double* STRESS,        // Stress
    double* STATEV,        // State variables
    double* DDSDDE,        // Material tangent
    double* SSE,           // Specific elastic energy
    double* SPD,           // Specific plastic dissipation
    double* SCD,           // Specific creep dissipation
    double* RPL,           // Volumetric heat generation per unit time 
    double* DDSDDT,        // Tangent (stress w.r.t temperature)
    double* DRPLDE,        // Tangent (RPL w.r.t. strain)
    double* DRPLDT,        // Tangent (RPL w.r.t. temperature)
    double* STRAN,         // Strain
    double* DSTRAN,        // Strain increment
    double* TIME,          // Time
    double* DTIME,         // Time increment
    double* TEMP,          // Temperature
    double* DTEMP,         // Temperature increment
    double* PREDEF,        // Pre-defined field variables
    double* DPRED,         // Pre-defined field variables increment
    char*   CMNAME,        // Material name
    int*    NDI,           // # of direct stress components
    int*    NSHR,          // # of shear stress components
    int*    NTENS,         // Size of stress/strain array
    int*    NSTATV,        // # of state variables
    double* PROPS,         // Material properties
    int*    NPROPS,        // # of material properties
    double* COORDS,        // Coordinates of this point
    double* DROT,          // (3,3) rotation increment matrix
    double* PNEWDT,        // Ratio of new to old step-size
    double* CELENT,        // Characteristic element length
    double* DFGRD0,        // (3,3) Deformation grad at start of increment
    double* DFGRD1,        // (3,3) Deformation grad at end of increment
    int*    NOEL,          // Element number
    int*    NPT,           // Integration point number
    int*    LAYER,         // Layer number
    int*    KSPT,          // Section point number within the current layer
    int*    KSTEP,         // Step number
    int*    KINC           // Increment number
    );

// Define the shared library as a macro
#define UMAT_SHARED_LIB "libumat.so"

// =======================================================================
//  class AbaqusMaterial
// =======================================================================

/** @brief 
 *  The class \c AbaqusMaterial implements an interface to Abaqus UMATs.
 *  In order to use this, first, the UMAT needs to be compiled as a 
 *  shared library (libumat.so). Check out this <a href="https://github.com/ritukeshbharali/Abaqus-UMAT-cpp" target="_blank">Github repo</a>
 *  to see how it is done.
 * 
 *  Below is an usage example, demonstrating how AbaqusMaterial is defined
 *  in the input file:
 * 
 *  \code
 *  material =
      {
        type     = "AbaqusUMAT";
        matProps = [100., 0.2];
        nStates  = 0;
        storeOldStrain = true;
        state  = "PLANE_STRAIN";  // not used
        rank   = 2;               // not used
      };
 *  \endcode
 * 
 *  Note that libumat.so must be presented in the same directory
 *  as the input file *.pro along with other input files.
 * 
 *  @warning Not thoroughly tested. Not all UMAT
 * function arguments are used.
 */


class AbaqusMaterial : public Material
{
 public:

  typedef  AbaqusMaterial Self;

  static const char*      UMAT_MATPROPS;
  static const char*      UMAT_NSTATES;
  static const char*      STORE_OLD_STRAIN;

  explicit                AbaqusMaterial

    ( int                   rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const;

  virtual void            updateConfig (); 

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const idx_t           ipoint );
  
  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const Vector&         dstrain,
      int                   ipoint );

  virtual void            allocPoints

    ( int                   ipCount );

  virtual Ref<Material>   clone  () const;

  virtual void            commit ();  
 
 protected:

  virtual                ~AbaqusMaterial   ();
  
 private:

  Vector                  matProps_; //!< Abaqus UMAT material props
  int                     nStates_;  //!< Abaqus UMAT number of states
  bool                    storeOld_; //!< Boolean flag to store old strain 
  umatFunc                umat;      //!< Function pointer to Abaqus UMAT
  IntVector               perm_;     //!< Abaqus to Jive stress/strain permutations

  // Store current and old step states (internal variables)

  class                   Hist_
  {
    public:

      Hist_( int n );

      Vector              v;
  };

  // Initialize state variables

  Flex<Hist_>             preState_;
  Flex<Hist_>             newState_;
  Flex<Hist_>*            latestState_;
  bool                    allocd_;

  // Initialize strains

  Flex<Hist_>             preStrain_;
  Flex<Hist_>             newStrain_;
  Flex<Hist_>*            latestStrain_;

};

#endif