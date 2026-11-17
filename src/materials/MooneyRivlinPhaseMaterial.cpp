/** @file MooneyRivlinPhaseMaterial.cpp
 *  @brief Compressible Mooney-Rivlin hyperelastic phase-fracture material.
 *
 *  Isochoric-volumetric strain energy:
 *
 *    Psi_dev = c1 (I1bar - 3) + c2 (I2bar - 3)
 *    Psi_vol = kappa/2 (ln J)^2
 *
 *  where
 *
 *    C     = I + 2 E
 *    I1    = tr(C)
 *    I2    = 1/2 [ I1^2 - tr(C^2) ]
 *    J     = sqrt(det(C))
 *    I1bar = J^(-2/3) I1
 *    I2bar = J^(-4/3) I2
 *
 *  The returned stress is second Piola-Kirchhoff stress and the
 *  returned tangent is dS/dE, for a Total-Lagrangian formulation.
 *
 *  Phase-field degradation:
 *
 *    S = gphi S_dev + S_vol
 *    D = gphi D_dev + D_vol
 */

#include <cmath>

#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/LUSolver.h>

#include "util/Constants.h"
#include "util/TensorUtils.h"
#include "MooneyRivlinPhaseMaterial.h"

using namespace jem;
using jem::numeric::LUSolver;

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* MooneyRivlinPhaseMaterial::C1_PROP    = "c1";
const char* MooneyRivlinPhaseMaterial::C2_PROP    = "c2";
const char* MooneyRivlinPhaseMaterial::KAPPA_PROP = "kappa";
const char* MooneyRivlinPhaseMaterial::STATE_PROP = "state";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

MooneyRivlinPhaseMaterial::MooneyRivlinPhaseMaterial

  ( const idx_t        rank,
    const Properties&  globdat )

  : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank == 2 || rank == 3 );

  c1_    = 0.5;
  c2_    = 0.0;
  kappa_ = 1.0;

  Psi_ = 0.0;

  derivPsi_.resize ( STRAIN_COUNTS[rank_] );
  derivPsi_ = 0.0;

  state_       = PlaneStrain;
  stateString_ = "PLANE_STRAIN";
}


MooneyRivlinPhaseMaterial::~MooneyRivlinPhaseMaterial ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void MooneyRivlinPhaseMaterial::configure

  ( const Properties&  props,
    const Properties&  globdat )
{
  props.get ( c1_,    C1_PROP    );
  props.get ( c2_,    C2_PROP    );
  props.get ( kappa_, KAPPA_PROP );

  if ( c1_ < 0.0 || c2_ < 0.0 )
  {
    throw Error (
      JEM_FUNC,
      "Mooney-Rivlin parameters c1 and c2 must be non-negative"
    );
  }

  if ( c1_ + c2_ <= 0.0 )
  {
    throw Error (
      JEM_FUNC,
      "Mooney-Rivlin requires c1 + c2 > 0"
    );
  }

  if ( kappa_ <= 0.0 )
  {
    throw Error (
      JEM_FUNC,
      "Mooney-Rivlin bulk modulus kappa must be positive"
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
        "Finite-strain Mooney-Rivlin plane stress is not implemented"
      );
    }
    else if ( stateString_ == "AXISYMMETRIC" )
    {
      throw Error (
        JEM_FUNC,
        "Axisymmetric Mooney-Rivlin formulation is not implemented"
      );
    }
    else
    {
      throw Error (
        JEM_FUNC,
        "Invalid state for MooneyRivlinPhaseMaterial: " + stateString_
      );
    }
  }
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void MooneyRivlinPhaseMaterial::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const
{
  conf.set ( C1_PROP,    c1_    );
  conf.set ( C2_PROP,    c2_    );
  conf.set ( KAPPA_PROP, kappa_ );

  if ( rank_ == 2 )
  {
    conf.set ( STATE_PROP, stateString_ );
  }
}

//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void MooneyRivlinPhaseMaterial::updateConfig ()
{
  // Nothing to update: c1_, c2_ and kappa_ are primary parameters.
}

//-----------------------------------------------------------------------
//   computeSplit_
//-----------------------------------------------------------------------

void MooneyRivlinPhaseMaterial::computeSplit_

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
  // 2D: {Exx,Eyy,Ezz,2Exy}
  // 3D: {Exx,Eyy,Ezz,2Exy,2Eyz,2Ezx}

  Matrix E = tensorUtils::voigt2tensorStrain ( strain );

  // Right Cauchy-Green tensor: C = I + 2E.

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
      "Mooney-Rivlin material encountered det(C) <= 0"
    );
  }

  const double J    = std::sqrt ( detC );
  const double lnJ  = std::log  ( J );
  const double I1   = tensorUtils::trace ( C );

  Matrix C2(3,3);
  C2 = 0.0;

  for ( idx_t i = 0; i < 3; ++i )
  {
    for ( idx_t j = 0; j < 3; ++j )
    {
      for ( idx_t k = 0; k < 3; ++k )
      {
        C2(i,j) += C(i,k) * C(k,j);
      }
    }
  }

  const double trC2 = tensorUtils::trace ( C2 );
  const double I2   = 0.5 * ( I1 * I1 - trC2 );

  const double Jm23 = std::pow ( J, -TWO_THIRD );
  const double Jm43 = Jm23 * Jm23;

  const double I1bar = Jm23 * I1;
  const double I2bar = Jm43 * I2;

  // ---------------------------------------------------------------
  // Isochoric/deviatoric PK2 stress.
  //
  // d(I1bar)/dC =
  //   J^(-2/3) [ I - I1/3 C^{-1} ]
  //
  // d(I2bar)/dC =
  //   J^(-4/3) [
  //     I1 I - C - 2 I2/3 C^{-1}
  //   ]
  //
  // S_dev = 2 c1 d(I1bar)/dC
  //       + 2 c2 d(I2bar)/dC
  // ---------------------------------------------------------------

  Matrix Sdev(3,3);
  Sdev = 0.0;

  for ( idx_t i = 0; i < 3; ++i )
  {
    for ( idx_t j = 0; j < 3; ++j )
    {
      const double dij = ( i == j ? 1.0 : 0.0 );

      const double A1 =
        Jm23 *
        (
          dij - ONE_THIRD * I1 * Cinv(i,j)
        );

      const double A2 =
        Jm43 *
        (
            I1 * dij
          - C(i,j)
          - TWO_THIRD * I2 * Cinv(i,j)
        );

      Sdev(i,j) = 2.0 * ( c1_ * A1 + c2_ * A2 );
    }
  }

  // ---------------------------------------------------------------
  // Volumetric PK2 stress.
  //
  // S_vol = kappa ln(J) C^{-1}
  // ---------------------------------------------------------------

  Matrix Svol(3,3);

  Svol = kappa_ * lnJ * Cinv;

  stressDev = tensorUtils::tensor2voigtStress ( Sdev, strCount );
  stressVol = tensorUtils::tensor2voigtStress ( Svol, strCount );

  // ---------------------------------------------------------------
  // Consistent isochoric/deviatoric tangent.
  //
  // Since S_dev = 2 d(Psi_dev)/dC and C = I + 2E,
  //
  //   D_dev = dS_dev/dE
  //         = 4 [ c1 H1 + c2 H2 ],
  //
  // where H1 and H2 are the Hessians of I1bar and I2bar
  // with respect to C.
  // ---------------------------------------------------------------

  Quadix Ddev(3,3,3,3);
  Ddev = 0.0;

  for ( idx_t i = 0; i < 3; ++i )
  {
    for ( idx_t j = 0; j < 3; ++j )
    {
      const double dij = ( i == j ? 1.0 : 0.0 );

      // Q_ij = I1 delta_ij - C_ij - 2/3 I2 Cinv_ij.
      const double Qij =
          I1 * dij
        - C(i,j)
        - TWO_THIRD * I2 * Cinv(i,j);

      for ( idx_t k = 0; k < 3; ++k )
      {
        for ( idx_t l = 0; l < 3; ++l )
        {
          const double dkl = ( k == l ? 1.0 : 0.0 );

          const double I4sym =
            0.5 *
            (
                ( i == k && j == l ? 1.0 : 0.0 )
              + ( i == l && j == k ? 1.0 : 0.0 )
            );

          const double CinvSym =
              Cinv(i,k) * Cinv(j,l)
            + Cinv(i,l) * Cinv(j,k);

          // Hessian of I1bar.

          const double H1 =
            Jm23 *
            (
              - ONE_THIRD *
                (
                    Cinv(k,l) * dij
                  + dkl * Cinv(i,j)
                )

              + ( I1 / 9.0 ) *
                  Cinv(i,j) * Cinv(k,l)

              + ( I1 / 6.0 ) * CinvSym
            );

          // Derivative of
          // Q_ij = I1 delta_ij - C_ij - 2/3 I2 Cinv_ij.

          const double dQ =
                dij * dkl
              - I4sym
              - TWO_THIRD *
                ( I1 * dkl - C(k,l) ) * Cinv(i,j)
              + ONE_THIRD * I2 * CinvSym;

          // Hessian of I2bar.

          const double H2 =
            Jm43 *
            (
              - TWO_THIRD * Cinv(k,l) * Qij
              + dQ
            );

          Ddev(i,j,k,l) =
            4.0 * ( c1_ * H1 + c2_ * H2 );
        }
      }
    }
  }

  // ---------------------------------------------------------------
  // Consistent volumetric tangent.
  //
  // D_vol_IJKL =
  //   kappa [
  //       Cinv_IJ Cinv_KL
  //     - ln(J) (Cinv_IK Cinv_JL + Cinv_IL Cinv_JK)
  //   ]
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

  // Only isochoric/deviatoric energy drives the phase field.

  Psi_ = c1_ * ( I1bar - 3.0 )
       + c2_ * ( I2bar - 3.0 );

  derivPsi_ = stressDev;
}

//-----------------------------------------------------------------------
//   update -- standard material response
//-----------------------------------------------------------------------

void MooneyRivlinPhaseMaterial::update

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

void MooneyRivlinPhaseMaterial::update

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

  stress = gphi * stressDev + stressVol;
  stiff  = gphi * stiffDev  + stiffVol;
}

//-----------------------------------------------------------------------
//   update -- split response
//-----------------------------------------------------------------------

void MooneyRivlinPhaseMaterial::update

  ( Vector&        stressP,
    Vector&        stressN,
    Matrix&        stiffP,
    Matrix&        stiffN,
    const Vector&  strain,
    int            ipoint )
{
  // P = degradable isochoric/deviatoric response.
  // N = non-degraded volumetric response.

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

double MooneyRivlinPhaseMaterial::givePsi () const
{
  return Psi_;
}

//-----------------------------------------------------------------------
//   giveDerivPsi
//-----------------------------------------------------------------------

Vector MooneyRivlinPhaseMaterial::giveDerivPsi () const
{
  return derivPsi_;
}

//-----------------------------------------------------------------------
//   clone
//-----------------------------------------------------------------------

Ref<Material> MooneyRivlinPhaseMaterial::clone () const
{
  return newInstance<MooneyRivlinPhaseMaterial> ( *this );
}
