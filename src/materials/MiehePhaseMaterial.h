
/** @file MiehePhaseMaterial.h
 *  @brief Implements a phase-field fracture material model with Miehe split.
 * 
 *  This class implements a phase-field fracture material
 *  model. Spectral decomposition of the strain energy 
 *  density is carried out based on Miehe et. al. (2010).
 *  (see DOI: 10.1016/j.cma.2010.04.011).
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 05 August 2022
 * 
 *  Updates (when, what and who)
 *     - [09 December 2022] Added Setters and a function updateConfig
 *       to update private members. (RB)
 */


#ifndef MIEHE_PHASEMATERIAL_H
#define MIEHE_PHASEMATERIAL_H

#include "util/TensorUtils.h"
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
//  class MiehePhaseMaterial
// =======================================================

// This class implements an isotropic phase-field fracture
// material

class MiehePhaseMaterial : public Material,
                           public PhaseFractureMaterial
{
 public:

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      STATE_PROP;

  explicit                MiehePhaseMaterial

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

  virtual                ~MiehePhaseMaterial   ();

 protected:

  // Elastic properties

  double                  area_;

 private:

  void                    getposProj_
  
    ( Quadix&             posProj,
      Matrix&             StrainT ); 

 private:

  ProblemType             state_;

  double                  K_;
  double                  G_;
  double                  lmda_;

  Vector                  stressP_;
  double                  Psi_;

  Matrix                  I2_;
  Matrix                  sigmaPos_,  sigmaNeg_;
  Matrix                  strainPos_, strainNeg_;

  Quadix                  posProj_, negProj_, Jacobian_, I4_;

};

//-----------------------------------------------------------------------
//   giveDim
//-----------------------------------------------------------------------

inline int    MiehePhaseMaterial::giveDim() const
{
  return rank_;
}

//-----------------------------------------------------------------------
//   giveStressP
//-----------------------------------------------------------------------

inline Vector  MiehePhaseMaterial::giveDerivPsi() const
{
  return stressP_;
}

//-----------------------------------------------------------------------
//   givePsi
//-----------------------------------------------------------------------

inline double  MiehePhaseMaterial::givePsi() const
{
  return Psi_;
}

//-----------------------------------------------------------------------
//   giveState
//-----------------------------------------------------------------------

inline ProblemType   MiehePhaseMaterial::giveState() const
{
  return state_;
}

#endif 
