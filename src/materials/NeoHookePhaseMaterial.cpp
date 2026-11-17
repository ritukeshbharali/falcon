/** @file NeoHookePhaseMaterial.cpp
 *  @brief Compressible Neo-Hookean hyperelastic phase-fracture material.
 *
 *  Isochoric-volumetric split:
 *
 *    Psi_dev = mu/2 ( J^(-2/3) I1 - 3 )
 *    Psi_vol = kappa/2 (ln J)^2
 *
 *    kappa = lambda + 2 mu / 3
 *
 *  with
 *
 *    C  = I + 2 E
 *    I1 = tr(C)
 *    J  = sqrt(det(C))
 *
 *  Second Piola-Kirchhoff stresses:
 *
 *    S_dev = mu J^(-2/3)
 *            [ I - (I1/3) C^(-1) ]
 *
 *    S_vol = kappa ln(J) C^(-1)
 *
 *  For phase-field degradation:
 *
 *    S = gphi S_dev + S_vol
 *    D = gphi D_dev + D_vol
 *
 *  Thus volumetric stiffness and stress are never degraded.
 */

#include <cmath>

#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/LUSolver.h>

#include "util/Constants.h"
#include "util/TensorUtils.h"
#include "NeoHookePhaseMaterial.h"

using namespace jem;
using jem::numeric::LUSolver;

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* NeoHookePhaseMaterial::LAMBDA_PROP = "lambda";
const char* NeoHookePhaseMaterial::MU_PROP     = "mu";
const char* NeoHookePhaseMaterial::STATE_PROP  = "state";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

NeoHookePhaseMaterial::NeoHookePhaseMaterial

  ( const idx_t        rank,
    const Properties&  globdat )

  : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank == 2 || rank == 3 );

  lambda_ = 1.0;
  mu_     = 1.0;
  kappa_  = lambda_ + TWO_THIRD * mu_;

  Psi_ = 0.0;

  derivPsi_.resize ( STRAIN_COUNTS[rank_] );
  derivPsi_ = 0.0;

  state_       = PlaneStrain;
  stateString_ = "PLANE_STRAIN";
}


NeoHookePhaseMaterial::~NeoHookePhaseMaterial ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void NeoHookePhaseMaterial::configure

  ( const Properties&  props,
    const Properties&  globdat )
{
  props.get ( lambda_, LAMBDA_PROP );
  props.get ( mu_,     MU_PROP     );

  if ( mu_ <= 0.0 )
  {
    throw Error ( JEM_FUNC, "mu must be positive" );
  }

  kappa_ = lambda_ + TWO_THIRD * mu_;

  if ( kappa_ <= 0.0 )
  {
    throw Error (
      JEM_FUNC,
      "bulk modulus kappa = lambda + 2/3 mu must be positive"
    );
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
        "Invalid state for NeoHookePhaseMaterial: " + stateString_
      );
    }
  }
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void NeoHookePhaseMaterial::getConfig

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

void NeoHookePhaseMaterial::updateConfig ()
{
  kappa_ = lambda_ + TWO_THIRD * mu_;
}

//-----------------------------------------------------------------------
//   computeSplit_
//-----------------------------------------------------------------------

void NeoHookePhaseMaterial::computeSplit_

  ( Vector&        stressDev,
    Vector&        stressVol,
    Matrix&        stiffDev,
    Matrix&        stiffVol,
    const Vector&  strain )
{
  const idx_t strCount = STRAIN_COUNTS[rank_];

  JEM_PRECHECK ( strain.size() == strCount );

  // Green-Lagrange strain tensor.
  //
  // 2D vector convention:
  //   { Exx, Eyy, Ezz, 2Exy }
  //
  // 3D vector convention:
  //   { Exx, Eyy, Ezz, 2Exy, 2Eyz, 2Ezx }

  Matrix E = tensorUtils::voigt2tensorStrain ( strain );

  // Right Cauchy-Green tensor:
  //
  // C = I + 2E

  Matrix C(3,3);

  C = 2.0 * E;

  for ( idx_t i = 0; i < 3; ++i )
  {
    C(i,i) += 1.0;
  }

  // C^{-1} and det(C).

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

  const double J     = std::sqrt ( detC );
  const double lnJ   = std::log  ( J );
  const double I1    = tensorUtils::trace ( C );
  const double Jm23  = std::pow ( J, -TWO_THIRD );
  const double I1d3  = ONE_THIRD * I1;

  // ---------------------------------------------------------------
  // Isochoric/deviatoric Second Piola-Kirchhoff stress
  //
  // S_dev = mu J^(-2/3) [ I - I1/3 C^{-1} ]
  // ---------------------------------------------------------------

  Matrix Sdev(3,3);

  Sdev = -mu_ * Jm23 * I1d3 * Cinv;

  for ( idx_t i = 0; i < 3; ++i )
  {
    Sdev(i,i) += mu_ * Jm23;
  }

  // ---------------------------------------------------------------
  // Volumetric Second Piola-Kirchhoff stress
  //
  // S_vol = kappa ln(J) C^{-1}
  // ---------------------------------------------------------------

  Matrix Svol(3,3);

  Svol = kappa_ * lnJ * Cinv;

  stressDev = tensorUtils::tensor2voigtStress ( Sdev, strCount );
  stressVol = tensorUtils::tensor2voigtStress ( Svol, strCount );

  // ---------------------------------------------------------------
  // Deviatoric consistent material tangent dS_dev / dE
  //
  // Ddev_IJKL = mu J^(-2/3) [
  //
  //    (I1/3) ( Cinv_IK Cinv_JL + Cinv_IL Cinv_JK )
  //
  //  + (2 I1/9) Cinv_IJ Cinv_KL
  //
  //  - (2/3) ( delta_IJ Cinv_KL + Cinv_IJ delta_KL )
  // ]
  // ---------------------------------------------------------------

  Quadix Ddev(3,3,3,3);
  Ddev = 0.0;

  for ( idx_t i = 0; i < 3; ++i )
  {
    for ( idx_t j = 0; j < 3; ++j )
    {
      const double dij = ( i == j ? 1.0 : 0.0 );

      for ( idx_t k = 0; k < 3; ++k )
      {
        for ( idx_t l = 0; l < 3; ++l )
        {
          const double dkl = ( k == l ? 1.0 : 0.0 );

          Ddev(i,j,k,l) =
            mu_ * Jm23 *
            (
                I1d3 *
                (
                    Cinv(i,k) * Cinv(j,l)
                  + Cinv(i,l) * Cinv(j,k)
                )

              + TWO_THIRD * I1d3 * Cinv(i,j) * Cinv(k,l)

              - TWO_THIRD *
                (
                    dij       * Cinv(k,l)
                  + Cinv(i,j) * dkl
                )
            );
        }
      }
    }
  }

  // ---------------------------------------------------------------
  // Volumetric consistent material tangent dS_vol / dE
  //
  // Dvol_IJKL = kappa [
  //
  //      Cinv_IJ Cinv_KL
  //
  //    - ln(J) ( Cinv_IK Cinv_JL + Cinv_IL Cinv_JK )
  // ]
  // ---------------------------------------------------------------

  Quadix Dvol(3,3,3,3);
  Dvol = 0.0;

  for ( idx_t i = 0; i < 3; ++i )
  {
    for ( idx_t j = 0; j < 3; ++j )
    {
      for ( idx_t k = 0; k < 3; ++k )
      {
        for ( idx_t l = 0; l < 3; ++l )
        {
          Dvol(i,j,k,l) =
            kappa_ *
            (
                Cinv(i,j) * Cinv(k,l)

              - lnJ *
                (
                    Cinv(i,k) * Cinv(j,l)
                  + Cinv(i,l) * Cinv(j,k)
                )
            );
        }
      }
    }
  }

  stiffDev = tensorUtils::tensor2voigt ( Ddev, strCount );
  stiffVol = tensorUtils::tensor2voigt ( Dvol, strCount );

  // Only the degradable part drives the phase field.

  Psi_ = 0.5 * mu_ * ( Jm23 * I1 - 3.0 );

  derivPsi_ = stressDev;
}

//-----------------------------------------------------------------------
//   update -- standard material response
//-----------------------------------------------------------------------

void NeoHookePhaseMaterial::update

  ( Vector&        stress,
    Matrix&        stiff,
    const Vector&  strain,
    const idx_t    ipoint )
{
  const idx_t strCount = STRAIN_COUNTS[rank_];

  Vector stressDev ( strCount );
  Vector stressVol ( strCount );

  Matrix stiffDev  ( strCount, strCount );
  Matrix stiffVol  ( strCount, strCount );

  computeSplit_ (
    stressDev,
    stressVol,
    stiffDev,
    stiffVol,
    strain
  );

  stress = stressDev + stressVol;
  stiff  = stiffDev  + stiffVol;
}

//-----------------------------------------------------------------------
//   update -- phase-field degraded response
//-----------------------------------------------------------------------

void NeoHookePhaseMaterial::update

  ( Vector&        stress,
    Matrix&        stiff,
    const Vector&  strain,
    const double&  gphi,
    int            ipoint )
{
  const idx_t strCount = STRAIN_COUNTS[rank_];

  Vector stressDev ( strCount );
  Vector stressVol ( strCount );

  Matrix stiffDev  ( strCount, strCount );
  Matrix stiffVol  ( strCount, strCount );

  computeSplit_ (
    stressDev,
    stressVol,
    stiffDev,
    stiffVol,
    strain
  );

  // Degradation acts only on the deviatoric/isochoric response.

  stress = gphi * stressDev + stressVol;
  stiff  = gphi * stiffDev  + stiffVol;
}

//-----------------------------------------------------------------------
//   update -- split response
//-----------------------------------------------------------------------

void NeoHookePhaseMaterial::update

  ( Vector&        stressP,
    Vector&        stressN,
    Matrix&        stiffP,
    Matrix&        stiffN,
    const Vector&  strain,
    int            ipoint )
{
  // For this material:
  //
  //   P = degradable deviatoric/isochoric contribution
  //   N = non-degraded volumetric contribution

  computeSplit_ (
    stressP,
    stressN,
    stiffP,
    stiffN,
    strain
  );
}

//-----------------------------------------------------------------------
//   givePsi
//-----------------------------------------------------------------------

double NeoHookePhaseMaterial::givePsi () const
{
  return Psi_;
}

//-----------------------------------------------------------------------
//   giveDerivPsi
//-----------------------------------------------------------------------

Vector NeoHookePhaseMaterial::giveDerivPsi () const
{
  return derivPsi_;
}

//-----------------------------------------------------------------------
//   clone
//-----------------------------------------------------------------------

Ref<Material> NeoHookePhaseMaterial::clone () const
{
  return newInstance<NeoHookePhaseMaterial> ( *this );
}
