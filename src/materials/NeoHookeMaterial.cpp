/** @file NeoHookeMaterial.cpp
 *  @brief Compressible Neo-Hookean hyperelastic material.
 *
 *  Strain measure:
 *      E = 1/2 (C - I)
 *
 *  Strain energy density:
 *      Psi = mu/2 (I1 - 3) - mu ln(J) + lambda/2 [ln(J)]^2
 *
 *  Second Piola-Kirchhoff stress:
 *      S = mu I + (lambda ln(J) - mu) C^{-1}
 *
 *  Consistent material tangent:
 *      C_IJKL = lambda Cinv_IJ Cinv_KL
 *             + (mu - lambda ln(J))
 *               (Cinv_IK Cinv_JL + Cinv_IL Cinv_JK)
 */

#include <cmath>

#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/LUSolver.h>

#include "util/Constants.h"
#include "util/TensorUtils.h"
#include "NeoHookeMaterial.h"

using namespace jem;
using jem::numeric::LUSolver;

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* NeoHookeMaterial::LAMBDA_PROP = "lambda";
const char* NeoHookeMaterial::MU_PROP     = "mu";
const char* NeoHookeMaterial::STATE_PROP  = "state";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

NeoHookeMaterial::NeoHookeMaterial

  ( const idx_t        rank,
    const Properties&  globdat )

  : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank == 2 || rank == 3 );

  lambda_ = 1.0;
  mu_     = 1.0;

  state_       = PlaneStrain;
  stateString_ = "PLANE_STRAIN";
}


NeoHookeMaterial::~NeoHookeMaterial ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void NeoHookeMaterial::configure

  ( const Properties&  props,
    const Properties&  globdat )
{
  props.get ( lambda_, LAMBDA_PROP );
  props.get ( mu_,     MU_PROP     );

  if ( mu_ <= 0.0 )
  {
    throw Error ( JEM_FUNC, "mu must be positive" );
  }

  if ( rank_ == 2 )
  {
    props.find ( stateString_, STATE_PROP );

    if ( stateString_ == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( stateString_ == "PLANE_STRESS" )
    {
      throw Error (
        JEM_FUNC,
        "Finite-strain Neo-Hookean plane stress is not implemented"
      );
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      throw Error (
        JEM_FUNC,
        "Axisymmetric Neo-Hookean formulation is not implemented"
      );
    }
    else
    {
      throw Error (
        JEM_FUNC,
        "Invalid state for NeoHookeMaterial: " + stateString_
      );
    }
  }
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void NeoHookeMaterial::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const
{
  conf.set ( LAMBDA_PROP, lambda_ );
  conf.set ( MU_PROP,     mu_     );

  if ( rank_ == 2 )
  {
    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void NeoHookeMaterial::updateConfig ()
{
  // Nothing to update. lambda_ and mu_ are primary material parameters.
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void NeoHookeMaterial::update

  ( Vector&        stress,
    Matrix&        stiff,
    const Vector&  strain,
    const idx_t    ipoint )
{
  const idx_t strCount = STRAIN_COUNTS[rank_];

  JEM_PRECHECK ( strain.size() == strCount );

  // Green-Lagrange strain:
  //
  // 2D: {Exx,Eyy,Ezz,2Exy}
  // 3D: {Exx,Eyy,Ezz,2Exy,2Eyz,2Ezx}
  //
  // The utility returns a full 3x3 tensor also for the 2D
  // {xx,yy,zz,xy} representation.

  Matrix E = tensorUtils::voigt2tensorStrain ( strain );

  // Right Cauchy-Green tensor: C = I + 2E

  Matrix C(3,3);

  C = 2.0 * E;

  for ( idx_t i = 0; i < 3; ++i )
  {
    C(i,i) += 1.0;
  }

  // Invert C. LUSolver::invert returns det(C) in detC.

  Matrix Cinv = C;

  double detC = 0.0;

  LUSolver::invert ( Cinv, detC );

  if ( detC <= 0.0 )
  {
    throw Error (
      JEM_FUNC,
      "Neo-Hookean material encountered det(C) <= 0"
    );
  }

  const double J   = std::sqrt ( detC );
  const double lnJ = std::log  ( J    );

  // Second Piola-Kirchhoff stress:
  //
  // S = mu I + (lambda ln(J) - mu) C^{-1}

  Matrix S(3,3);

  S = ( lambda_ * lnJ - mu_ ) * Cinv;

  for ( idx_t i = 0; i < 3; ++i )
  {
    S(i,i) += mu_;
  }

  stress = tensorUtils::tensor2voigtStress ( S, strCount );

  // Consistent material tangent:
  //
  // dS_IJ/dE_KL =
  // lambda Cinv_IJ Cinv_KL
  // + (mu - lambda lnJ)
  //   (Cinv_IK Cinv_JL + Cinv_IL Cinv_JK)

  Quadix Ct(3,3,3,3);

  Ct = 0.0;

  const double a = mu_ - lambda_ * lnJ;

  for ( idx_t i = 0; i < 3; ++i )
  {
    for ( idx_t j = 0; j < 3; ++j )
    {
      for ( idx_t k = 0; k < 3; ++k )
      {
        for ( idx_t l = 0; l < 3; ++l )
        {
          Ct(i,j,k,l) =
              lambda_ * Cinv(i,j) * Cinv(k,l)
            + a * (
                  Cinv(i,k) * Cinv(j,l)
                + Cinv(i,l) * Cinv(j,k)
              );
        }
      }
    }
  }

  stiff = tensorUtils::tensor2voigt ( Ct, strCount );
}

//-----------------------------------------------------------------------
//   clone
//-----------------------------------------------------------------------

Ref<Material> NeoHookeMaterial::clone () const
{
  return newInstance<NeoHookeMaterial> ( *this );
}
