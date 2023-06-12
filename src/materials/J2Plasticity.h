
/** @file J2Plasticity.h
 *  @brief Implements a mixed isotropic-kinemetic linear hardening J2
 *  plasticity model.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 12 June 2023
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2023] 
 */

#ifndef J2_PLASTICTY_H
#define J2_PLASTICTY_H

#include <jem/util/Flex.h>

#include "util/BasicUtils.h"
#include "util/Arrays.h"
#include "Material.h"
#include "HookeMaterial.h"

using jem::util::Flex;

/*
enum ProblemType {
  PlaneStrain,
  PlaneStress,
  AxiSymmetric
};
*/

// =======================================================
//  class J2Plasticity
// =======================================================

class J2Plasticity : public Material
{
 public:

  typedef J2Plasticity    Self;
  typedef Material        Super;

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      YIELD_STRESS_PROP;
  static const char*      H_MODULUS_PROP;
  static const char*      R_PROP;
  static const char*      STATE_PROP;

  explicit                J2Plasticity

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

  inline void            getHistoryNames   
    
    ( const StringVector&   hnames ) const;

 protected:

  virtual                ~J2Plasticity   ();

 protected:

  // Hooke material

  Ref<HookeMaterial>      elasticMat_;
  Matrix                  elasticStiffMat_;

  // History

  class                   Hist_
  {
    public:

      Hist_( );

      Vector              sig;       // stress
      Vector              epsp;      // plastic strain
      Vector              alf;       // back stress
      double              kap;       // hardening variable
      bool                loading;   // plastic loading  
  };

  Flex<Hist_>             preHist_;
  Flex<Hist_>             newHist_;
  Flex<Hist_>*            latestHist_;

 private:

  double                  sigY_;
  double                  H_;
  double                  r_;

  // Derived

  double                  K_;
  double                  G_;

  Matrix                  P_;
  Matrix                  Q_;

};

#endif 
