
/** @file AmorPhaseMaterial.h
 *  @brief Implements a phase-field fracture material model with Amor split.
 *  
 *  This class implements a phase-field fracture material
 *  model. Vol-dev decomposition of the strain energy 
 *  density is carried out based on Amor et. al. (2009).
 *  (see DOI: 10.1016/j.jmps.2009.04.011).
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 05 August 2022
 * 
 *  Updates (when, what and who)
 *     - [09 December 2022] Added Setters and a function updateConfig
 *       to update private members. (RB)
 */

#ifndef AMOR_PHASEMATERIAL_H
#define AMOR_PHASEMATERIAL_H

#include "Material.h"
#include "PhaseFractureMaterial.h"

/*
enum ProblemType {
  PlaneStrain,
  PlaneStress,
  AxiSymmetric
};
*/

// =======================================================
//  class AmorPhaseMaterial
// =======================================================

class AmorPhaseMaterial : public Material,
                          public PhaseFractureMaterial
{
 public:

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      STATE_PROP;

  explicit                AmorPhaseMaterial

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
      const double&         gphi,
      int                   ipoint );

  virtual void            update

    ( Vector&               stressP,
      Vector&               stressN,
      Matrix&               stiffP,
      Matrix&               stiffN,
      const Vector&         strain,
      int                   ipoint );

  virtual Ref<Material>   clone ( ) const;         

  //  Getters (provides read-only/copy access private members of the class)  

  inline int              giveDim           () const;
  inline ProblemType      giveState         () const;
  inline Vector           giveDerivPsi      () const;
  inline double           givePsi           () const;

 protected:

  virtual                ~AmorPhaseMaterial   ();

 protected:

  // Elastic properties

  double                  K_;
  double                  G_;
  double                  area_;

 private:

  ProblemType             state_;
  Vector                  stressP_;
  double                  Psi_;

  // Projection operators
  
  Matrix  PVol_;
  Matrix  PDev_;

};

//-----------------------------------------------------------------------
//   giveDim
//-----------------------------------------------------------------------

inline int    AmorPhaseMaterial::giveDim() const
{
  return rank_;
}

//-----------------------------------------------------------------------
//   giveDerivPsi
//-----------------------------------------------------------------------

inline Vector  AmorPhaseMaterial::giveDerivPsi() const
{
  return stressP_;
}

//-----------------------------------------------------------------------
//   givePsi
//-----------------------------------------------------------------------

inline double  AmorPhaseMaterial::givePsi() const
{
  return Psi_;
}

//-----------------------------------------------------------------------
//   giveState
//-----------------------------------------------------------------------

inline ProblemType   AmorPhaseMaterial::giveState() const
{
  return state_;
}

#endif 
