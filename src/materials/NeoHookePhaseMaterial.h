/** @file NeoHookePhaseMaterial.h
 *  @brief Compressible Neo-Hookean hyperelastic phase-fracture material.
 *
 *  Uses an isochoric-volumetric split:
 *
 *    Psi = Psi_dev + Psi_vol
 *
 *    Psi_dev = mu/2 ( J^(-2/3) I1 - 3 )
 *    Psi_vol = kappa/2 (ln J)^2
 *
 *  with kappa = lambda + 2 mu / 3.
 *
 *  In the phase-field update only the deviatoric/isochoric contribution
 *  is degraded:
 *
 *    S = gphi S_dev + S_vol
 *    C = gphi C_dev + C_vol
 */

#ifndef NEOHOOKEPHASE_MATERIAL_H
#define NEOHOOKEPHASE_MATERIAL_H

#include "Material.h"
#include "PhaseFractureMaterial.h"

class NeoHookePhaseMaterial : public Material,
                         public PhaseFractureMaterial
{
 public:

  static const char*      LAMBDA_PROP;
  static const char*      MU_PROP;
  static const char*      STATE_PROP;

  explicit                NeoHookePhaseMaterial

    ( const idx_t          rank,
      const Properties&    globdat );

  virtual void            configure

    ( const Properties&    props,
      const Properties&    globdat );

  virtual void            getConfig

    ( const Properties&    conf,
      const Properties&    globdat ) const;

  virtual void            updateConfig ();

  // Standard material update: undegraded Neo-Hookean response.

  virtual void            update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint );

  // Phase-field update: degradation acts only on deviatoric response.

  virtual void            update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const double&        gphi,
      int                  ipoint );

  // Split response for micromorphic phase-field formulations:
  // stressP/stiffP = degradable deviatoric part
  // stressN/stiffN = non-degraded volumetric part.

  virtual void            update

    ( Vector&              stressP,
      Vector&              stressN,
      Matrix&              stiffP,
      Matrix&              stiffN,
      const Vector&        strain,
      int                  ipoint );

  virtual double          givePsi      () const;
  virtual Vector          giveDerivPsi () const;

  virtual Ref<Material>   clone () const;

 protected:

  virtual                ~NeoHookePhaseMaterial ();

 private:

  void                    computeSplit_

    ( Vector&              stressDev,
      Vector&              stressVol,
      Matrix&              stiffDev,
      Matrix&              stiffVol,
      const Vector&        strain );

 private:

  double                  lambda_;
  double                  mu_;
  double                  kappa_;

  double                  Psi_;
  Vector                  derivPsi_;

  String                  stateString_;
  ProblemType             state_;
};

#endif
