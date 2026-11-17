/** @file NeoHookeMaterial.h
 *  @brief Compressible Neo-Hookean hyperelastic material.
 *
 *  The material update takes Green-Lagrange strain and returns
 *  second Piola-Kirchhoff stress and the consistent material tangent
 *  dS/dE for use in a Total-Lagrangian finite-strain formulation.
 */

#ifndef NEOHOOKE_MATERIAL_H
#define NEOHOOKE_MATERIAL_H

#include "Material.h"

class NeoHookeMaterial : public Material
{
 public:

  static const char*      LAMBDA_PROP;
  static const char*      MU_PROP;
  static const char*      STATE_PROP;

  explicit                NeoHookeMaterial

    ( const idx_t          rank,
      const Properties&    globdat );

  virtual void            configure

    ( const Properties&    props,
      const Properties&    globdat );

  virtual void            getConfig

    ( const Properties&    conf,
      const Properties&    globdat ) const;

  virtual void            updateConfig ();

  virtual void            update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint );

  virtual Ref<Material>    clone () const;

 protected:

  virtual                ~NeoHookeMaterial ();

 private:

  double                  lambda_;
  double                  mu_;

  String                  stateString_;
  ProblemType             state_;
};

#endif
