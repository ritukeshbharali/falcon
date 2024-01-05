
/** @file MiehePhaseMaterial.cpp
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

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/EigenUtils.h>

#include "util/Constants.h"
#include "util/MathUtils.h"
#include "util/TensorUtils.h"

#include "MiehePhaseMaterial.h"

using namespace jem;
using namespace jem::io;
using namespace tensorUtils;

using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  MiehePhaseMaterial::YOUNG_PROP          = "young";
const char*  MiehePhaseMaterial::POISSON_PROP        = "poisson";
const char*  MiehePhaseMaterial::STATE_PROP          = "state";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


MiehePhaseMaterial::MiehePhaseMaterial 

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  young_   = 1.0;
  poisson_ = 1.0;

  stressP_.resize( STRAIN_COUNTS[rank_] );
  stressP_ = 0.0;
}


MiehePhaseMaterial::~MiehePhaseMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void MiehePhaseMaterial::configure ( const Properties& props,
                                const Properties& globdat )
{
  using jem::maxOf;

  props.find ( young_,   YOUNG_PROP,   0.0, maxOf( young_ ) );
  props.find ( poisson_, POISSON_PROP, 0.0, 0.5 );
  
  double K, G;

  if ( ( props.find ( K, "bulk_modulus"  ) ) && 
       ( props.find ( G, "shear_modulus" ) ) )
  {
    young_   = 9*K*G / ( 3*K + G );
    poisson_ = ( 3*K - 2*G ) / 2. / ( 3*K + G );
  }

  K_    = young_ / ( 3.0 * ( 1.0 - 2.0 * poisson_ ) );
  G_    = young_ / ( 2.0 * ( 1.0 + poisson_) );
  lmda_ = K_ - 0.666666667 * G_;

  if ( rank_ == 1 )
  {
    props.get ( area_, "area" );
  }

  // read problem type, plane stress ...

  String state;

  if ( rank_ == 2  )
  {
    props.get( state, STATE_PROP );

    if      ( state == "PLANE_STRAIN" )
    {
      state_ = PlaneStrain;
    }
    else if ( state == "PLANE_STRESS" )
    {
      state_ = PlaneStress;
    }
    else if ( state == "AXISYMMETRIC" )
    {
      state_ = AxiSymmetric;
    }
  }

  // Initialize tensors

  I2_.resize(3,3);
  I4_.resize(3,3,3,3);

  I2_ = 0.0;
  for ( int i = 0; i < 3; ++i )
  {
    I2_(i,i) = 1.0;
  }

  I4_ = getSymI4 ( 3 );

  sigmaPos_.resize(3,3);
  sigmaNeg_.resize(3,3);
  strainPos_.resize(3,3);
  strainNeg_.resize(3,3);

  posProj_ .resize(3,3,3,3);
  negProj_ .resize(3,3,3,3);
  Jacobian_.resize(3,3,3,3);
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void MiehePhaseMaterial::getConfig ( const Properties& conf, 
                                const Properties& globdat ) const
{
  conf.set ( YOUNG_PROP,   young_   );
  conf.set ( POISSON_PROP, poisson_ );
  conf.set ( STATE_PROP,   state_   );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void MiehePhaseMaterial::updateConfig ()

{
  K_    = young_ / ( 3.0 * ( 1.0 - 2.0 * poisson_ ) );
  G_    = young_ / ( 2.0 * ( 1.0 + poisson_) );
  lmda_ = K_ - 0.666666667 * G_;
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void MiehePhaseMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{
  // do nothing ... this should never be called from a phase-field
  // FE model.
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void MiehePhaseMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const double&         gphi,
      int                   ipoint )

{

  // Fill 3D strain ( voigt to tensor)

  Matrix strainT = voigt2tensorStrain( strain );

  // Compute the positve and negative projection tensors

  getposProj_( posProj_, strainT );
  negProj_ = I4_ - posProj_;

  // Compute positive and negative strains

  strainPos_ = doubleDot ( posProj_, strainT );
  strainNeg_ = strainT - strainPos_;

  double straintrace = trace( strainT );

  // Compute the positive and negative stress tensors

  sigmaPos_ = lmda_ *  macaulayP ( straintrace ) * I2_ + 2.0 * G_ * strainPos_;
  sigmaNeg_ = lmda_ *  macaulayN ( straintrace ) * I2_ + 2.0 * G_ * strainNeg_;

  // Transform stress tensors to voigt vectors

  Vector stressP = tensor2voigtStress( sigmaPos_, STRAIN_COUNTS[rank_] );
  Vector stressN = tensor2voigtStress( sigmaNeg_, STRAIN_COUNTS[rank_] );

  stress = gphi * stressP + stressN;

  // Store current positive stress and energy

  stressP_ = stressP;

  Psi_ = 0.5 * lmda_ * ::pow( macaulayP( straintrace ), 2.0  ) 
                   + G_ * doubleDot( strainPos_, strainPos_ );

  // Compute the material tangent stiffness matrix

  double signpos = heavisideP ( straintrace );
  double signneg = 1.0 - signpos;

  Quadix jac (3,3,3,3);
  jac = 0.0;

  jac = gphi * ( 2.0 * G_ * posProj_ + lmda_ * signpos * otimes( I2_, I2_ )  )
             + ( 2.0 * G_ * negProj_ + lmda_ * signneg * otimes( I2_, I2_ )  );

  stiff = tensor2voigt( jac , STRAIN_COUNTS[rank_] );

}

//-----------------------------------------------------------------------
//   update (overloaded)
//-----------------------------------------------------------------------

void MiehePhaseMaterial::update

    ( Vector&               stressP,
      Vector&               stressN,
      Matrix&               stiffP,
      Matrix&               stiffN,
      const Vector&         strain,
      int                   ipoint )

{

  // Fill 3D strain ( voigt to tensor)

  Matrix strainT = voigt2tensorStrain( strain );

  // Compute the positve and negative projection tensors

  getposProj_( posProj_, strainT );
  negProj_ = I4_ - posProj_;

  // Compute positive and negative strains

  strainPos_ = doubleDot ( posProj_, strainT );
  strainNeg_ = strainT - strainPos_;

  double straintrace = trace( strainT );

  // Compute the positive and negative stress tensors

  sigmaPos_ = lmda_ *  macaulayP ( straintrace ) * I2_ + 2.0 * G_ * strainPos_;
  sigmaNeg_ = lmda_ *  macaulayN ( straintrace ) * I2_ + 2.0 * G_ * strainNeg_;

  // Transform stress tensors to voigt vectors

  stressP = tensor2voigtStress( sigmaPos_, STRAIN_COUNTS[rank_] );
  stressN = tensor2voigtStress( sigmaNeg_, STRAIN_COUNTS[rank_] );

  // Store current positive stress and energy

  stressP_ = stressP;

  Psi_ = 0.5 * lmda_ * ::pow( macaulayP( straintrace ), 2.0  ) 
                   + G_ * doubleDot( strainPos_, strainPos_ );

  // Compute the material tangent stiffness matrix

  double signpos = heavisideP ( straintrace );
  double signneg = 1.0 - signpos;

  Quadix jacP (3,3,3,3);
  Quadix jacN (3,3,3,3);

  jacP = jacN = 0.0;

  jacP = ( 2.0 * G_ * posProj_ + lmda_ * signpos * otimes( I2_, I2_ )  );
  jacN = ( 2.0 * G_ * negProj_ + lmda_ * signneg * otimes( I2_, I2_ )  );

  stiffP = tensor2voigt( jacP , STRAIN_COUNTS[rank_] );
  stiffN = tensor2voigt( jacN , STRAIN_COUNTS[rank_] );

}


//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> MiehePhaseMaterial::clone ( ) const

{
  return newInstance<MiehePhaseMaterial> ( *this );
}


//-----------------------------------------------------------------------
//   getposProj_
//-----------------------------------------------------------------------

void MiehePhaseMaterial::getposProj_

  ( Quadix&  posProj,
    Matrix&  StrainT )

{
  Matrix  eigVecs;    eigVecs .resize ( 3, 3 );
  Vector  eigVals;    eigVals .resize ( 3 );
  Vector  epos;       epos    .resize ( 3 );
  Vector  diag;       diag    .resize ( 3 );
  Matrix  Ma;         Ma      .resize ( 3, 3 );
  Matrix  Mb;         Mb      .resize ( 3, 3 );
  Quadix  Gab;        Gab     .resize ( 3 , 3, 3, 3 );
  Quadix  Gba;        Gba     .resize ( 3 , 3, 3, 3 );
  double  theta_ab;
  
  posProj  = 0.0;
  Ma       = 0.0;
  Mb       = 0.0;

  // Compute eigen values and vectors of B

  using jem::numeric::EigenUtils;

  EigenUtils::symSolve ( eigVals, eigVecs, StrainT );

  for ( int i = 0; i < 3 ; i++ )
  {
     epos[i] = macaulayP( eigVals[i] );
     diag[i] = 0.0;
     if ( eigVals[i]>0.0 )
     {
        diag[i]=1.0;
     }
  }
  
  for ( int iRnk = 0; iRnk < 3 ; iRnk++ )
  {
     Ma       = otimes(eigVecs(ALL,iRnk),eigVecs(ALL,iRnk));
     posProj += otimes(Ma,Ma)*diag[iRnk];  
  }

  for ( int a = 0; a < 3 ; a++ )
  {
    for ( int b =0; b < a ; b++ )
    {
       Ma = otimes(eigVecs(ALL,a),eigVecs(ALL,a));
       Mb = otimes(eigVecs(ALL,b),eigVecs(ALL,b));
       
       Gab = otimesu(Ma,Mb) + otimesl(Ma,Mb);
       Gba = otimesu(Mb,Ma) + otimesl(Mb,Ma);
       
       if ( fabs( eigVals[a] - eigVals[b] ) <= 1e-13)
       {
         theta_ab = 0.5*(diag[a]+diag[b])/2.0;
       }
       else
       {
         theta_ab = 0.5*(epos[a]-epos[b])/(eigVals[a]-eigVals[b]);
       }
       posProj += theta_ab*(Gab+Gba);
    }
  }
  
}