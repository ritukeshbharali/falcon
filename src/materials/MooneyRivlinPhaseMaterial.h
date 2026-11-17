/** @file MooneyRivlinPhaseMaterial.h
 *  @brief Compressible Mooney-Rivlin hyperelastic phase-fracture material.
 *
 *  Uses an isochoric-volumetric split:
 *
 *    Psi = Psi_dev + Psi_vol
 *
 *    Psi_dev = c1 (I1bar - 3) + c2 (I2bar - 3)
 *    Psi_vol = kappa/2 (ln J)^2
 *
 *  with
 *
 *    I1bar = J^(-2/3) I1
 *    I2bar = J^(-4/3) I2.
 *
 *  In the phase-field update only the isochoric/deviatoric response
 *  is degraded:
 *
 *    S = gphi S_dev + S_vol
 *    D = gphi D_dev + D_vol
 */

#ifndef MOONEYRIVLINPHASE_MATERIAL_H
#define MOONEYRIVLINPHASE_MATERIAL_H

#include "Material.h"
#include "PhaseFractureMaterial.h"

class MooneyRivlinPhaseMaterial : public Material,
                             public PhaseFractureMaterial
{
 public:

  static const char*      C1_PROP;
  static const char*      C2_PROP;
  static const char*      KAPPA_PROP;
  static const char*      STATE_PROP;

  explicit                MooneyRivlinPhaseMaterial

    ( const idx_t          rank,
      const Properties&    globdat );

  virtual void            configure

    ( const Properties&    props,
      const Properties&    globdat );

  virtual void            getConfig

    ( const Properties&    conf,
      const Properties&    globdat ) const;

  virtual void            updateConfig ();

  // Standard undegraded response.

  virtual void            update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint );

  // Phase-field response: only isochoric/deviatoric part degraded.

  virtual void            update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const double&        gphi,
      int                  ipoint );

  // Split response:
  // P = degradable isochoric/deviatoric contribution
  // N = non-degraded volumetric contribution.

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

  virtual                ~MooneyRivlinPhaseMaterial ();

 private:

  void                    computeSplit_

    ( Vector&              stressDev,
      Vector&              stressVol,
      Matrix&              stiffDev,
      Matrix&              stiffVol,
      const Vector&        strain );

 private:

  double                  c1_;
  double                  c2_;
  double                  kappa_;

  double                  Psi_;
  Vector                  derivPsi_;

  String                  stateString_;
  ProblemType             state_;
};

#endif
