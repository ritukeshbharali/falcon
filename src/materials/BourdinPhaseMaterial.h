
/** @file BourdinPhaseMaterial.h
 *  @brief Implements a phase-field fracture material model with no split.
 * 
 *  This class implements a phase-field fracture material
 *  model without any split in the strain energy density
 *  (see DOI: 10.1016/S0022-5096(99)00028-9).
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 05 August 2022
 * 
 *  Updates (when, what and who)
 *     - [09 December 2022] Added Setters and a function updateConfig
 *       to update private members. (RB)
 */

#ifndef BOURDIN_PHASEMATERIAL_H
#define BOURDIN_PHASEMATERIAL_H

#include "util/BasicUtils.h"
#include "Material.h"
#include "HookeMaterial.h"
#include "PhaseFractureMaterial.h"

/*
enum ProblemType {
  PlaneStrain,
  PlaneStress,
  AxiSymmetric
};
*/

// =======================================================
//  class BourdinPhaseMaterial
// =======================================================

class BourdinPhaseMaterial : public Material,
                             public PhaseFractureMaterial
{
 public:

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      STATE_PROP;

  explicit                BourdinPhaseMaterial

    ( int                   rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const;

  virtual void            updateConfig ();    

  virtual void           update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint );

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const double&         gphi,
      int                   ipoint );

  virtual void            update

    ( Vector&               stressP,
      Vector&               stressN,
      Matrix&               stiffP,
      Matrix&               stiffN,
      const Vector&         strain,
      int                   ip     ); 

  virtual Ref<Material>   clone ( ) const;     

  inline Vector           giveDerivPsi      () const;
  inline double           givePsi           () const;

 protected:

  virtual                ~BourdinPhaseMaterial   ();

 protected:

  // Hooke material

  Ref<HookeMaterial>      elasticMat_;
  Matrix                  elasticStiffMat_;

 private:

  Vector                  stressP_;
  double                  Psi_;

};

//-----------------------------------------------------------------------
//   giveStressP
//-----------------------------------------------------------------------

inline Vector  BourdinPhaseMaterial::giveDerivPsi() const
{
  return stressP_;
}

//-----------------------------------------------------------------------
//   givePsi
//-----------------------------------------------------------------------

inline double  BourdinPhaseMaterial::givePsi() const
{
  return Psi_;
}

#endif 
