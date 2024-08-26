
/** @file J2Material.h
 *  @brief Implements a mixed isotropic-kinemetic linear 
 *  hardening J2 plasticity model.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 19 Aug 2024
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2024] 
 */

#ifndef J2MATERIAL_H
#define J2MATERIAL_H

#include <jem/util/Flex.h>

#include "Material.h"

using jem::String;
using jem::util::Flex;

/*
enum ProblemType {
  PlaneStrain,
  PlaneStress,
  AxiSymmetric
};
*/

// =======================================================
//  class J2Material
// =======================================================

/** @brief Implements a J2 (von Mises) material model.
 * 
 *  The class \c J2Material implements a J2 (von Mises)
 *  plasticity material model, with mixed isotropic and
 *  kinematic hardening. The hardening is linear, as such
 *  the return to the yield surface is explicitly computed.
 * 
 *  Below is an example how the material model is defined:
 * 
 *  \code
 *  material =
    {
      type   =  "J2";            // Material type
      rank   =  2;               // Rank (2D in this case!)
      state  =  "PLANE_STRAIN";  // Material state

      young     = 100.e+09;      // Young's modulus
      poisson   = 0.3;           // Poisson ratio
      rho       = 0.0;           // Material density
      
      sig0      = 1.0;           // Initial yield stress
      h         = 1.0;           // Hardening modulus
      r         = 0.5;           // Interaction parameter
    };
 *  \endcode 
 */

class J2Material : public Material
{
 public:

  typedef J2Material      Self;
  typedef Material        Super;

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      YIELD_STRESS0_PROP;
  static const char*      HARD_MODULUS_PROP;
  static const char*      R_PROP;
  static const char*      STATE_PROP;

  explicit                J2Material

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

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint );

  virtual void           allocPoints

    ( const idx_t          npoints );  

  virtual Ref<Material>   clone ( ) const;     

  virtual void            commit ();

  //  Getters (provides read-only/copy access private/protected members of the class)

  virtual void           getHistory

    ( Vector&              hvals,
      const idx_t          mpoint );

 protected:

  virtual                ~J2Material   ();

 protected:

  // History

  class                   Hist_
  {
    public:

      Hist_( );
      void                toVector ( const Vector& v ) const;
      void                toPrint  () const;

      Vector              sig;         //!< Stress
      Vector              eps;         //!< Strain
      Vector              alf;         //!< Back stress
      double              kap;         //!< Hardening variable
      double              mu;          //!< Integrated plastic multiplier
      bool                loading;     //!< Plastic loading
  };

  Flex<Hist_>             preHist_;    //!< Prev. step history
  Flex<Hist_>             newHist_;    //!< Curr. step history
  Flex<Hist_>*            latestHist_; //!< Ptr. to latest history

 private:

  // Input parameters
  double                  sigY_;       //!< Initial yield stress
  double                  h_;          //!< Hardening modulus
  double                  r_;          //!< Interaction parameter

  // Derived parameters  
  double                  strCount_;   //!< Strain count
  double                  K_;          //!< Bulk modulus
  double                  G_;          //!< Shear modulus
  
  // Projection matrices
  Matrix                  PVol_;       //!< Strain volumetric projection
  Matrix                  PDev_;       //!< Strain deviatoric projection
  Matrix                  QDev_;       //!< Stress deviatoric projection   
  Matrix                  stiff0_;     //!< Elastic stiffness
};

#endif