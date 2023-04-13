
/** @file AmorPhaseMaterial.cpp
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


#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/utilities.h>

#include "util/BasicUtils.h"
#include "AmorPhaseMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  AmorPhaseMaterial::YOUNG_PROP          = "young";
const char*  AmorPhaseMaterial::POISSON_PROP        = "poisson";
const char*  AmorPhaseMaterial::STATE_PROP          = "state";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


AmorPhaseMaterial::AmorPhaseMaterial 

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  young_   = 1.0;
  poisson_ = 1.0;

  stressP_.resize( STRAIN_COUNTS[rank_] );
  stressP_ = 0.0;
}


AmorPhaseMaterial::~AmorPhaseMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void AmorPhaseMaterial::configure ( const Properties& props,
                                const Properties& globdat )
{
  using jem::maxOf;

  props.find ( young_,   YOUNG_PROP,   0.0, maxOf( young_ ) );
  props.find ( poisson_, POISSON_PROP, 0.0, 0.5 );

  K_ = young_ / ( 3.0 * ( 1.0 - 2.0 * poisson_ ) );
  G_ = young_ / ( 2.0 * ( 1.0 + poisson_) );
  
  if ( ( props.find ( K_, "bulk_modulus"  ) ) && 
       ( props.find ( G_, "shear_modulus" ) ) )
  {
    young_   = 9*K_*G_ / ( 3*K_ + G_ );
    poisson_ = ( 3*K_ - 2*G_ ) / 2. / ( 3*K_ + G_ );
  }

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

  // compute volumetric and deviatoric projection matrices

  PVol_.resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  PVol_ = 0.0;

  PDev_.resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  PDev_ = 0.0;

  for ( int i = 0; i<3; i++)
  {
    for ( int j = 0; j<3; j++)
    {
      PVol_(i,j) = 1.0;
      i == j ? PDev_(i,j) = 2.0/3.0: PDev_(i,j) = -1.0/3.0;
    }
  }

  if ( rank_ == 2 )
  {
    PDev_(3,3) = 0.5;
  }
  else if ( rank_ == 3 )
  {
    PDev_(3,3) = 0.5;
    PDev_(4,4) = 0.5;
    PDev_(5,5) = 0.5;
  }

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void AmorPhaseMaterial::getConfig ( const Properties& conf, 
                                    const Properties& globdat ) const
{
  conf.set ( YOUNG_PROP,   young_   );
  conf.set ( POISSON_PROP, poisson_ );
  conf.set ( STATE_PROP,   state_   );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void AmorPhaseMaterial::updateConfig ()

{
  K_ = young_ / ( 3.0 * ( 1.0 - 2.0 * poisson_ ) );
  G_ = young_ / ( 2.0 * ( 1.0 + poisson_) );
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void AmorPhaseMaterial::update

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

void AmorPhaseMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const double&         gphi,
      int                   ipoint )

{

  JEM_ASSERT( strain.size() > 3 );

  double traceStrain = strain[0] + strain[1] + strain[2];

  stiff = 2.0 * G_ * PDev_;

  if ( traceStrain > 0.0 )
  {
    stiff   += K_ * PVol_;
  }

  matmul ( stressP_, stiff, strain );

  Psi_   = 0.5 * dot ( stressP_,strain );
  stiff *= gphi;

  if ( traceStrain <= 0.0 )
  {
    stiff   += K_ * PVol_;
  }

  matmul ( stress, stiff, strain );
}


//-----------------------------------------------------------------------
//   update (overloaded)
//-----------------------------------------------------------------------

void AmorPhaseMaterial::update

    ( Vector&               stressP,
      Vector&               stressN,
      Matrix&               stiffP,
      Matrix&               stiffN,
      const Vector&         strain,
      int                   ipoint )

{

  JEM_ASSERT( strain.size() > 3 );

  double traceStrain = strain[0] + strain[1] + strain[2];

  stiffP = 2.0 * G_ * PDev_;

  if ( traceStrain > 0.0 )
  {
    stiffP   += K_ * PVol_;
  }
  else
  {
    stiffN    = K_ * PVol_;
  }

  matmul ( stressP, stiffP, strain );
  matmul ( stressN, stiffN, strain );

  stressP_ = stressP; 
  Psi_     = 0.5 * dot ( stressP_,strain );

}


//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> AmorPhaseMaterial::clone ( ) const

{
  return newInstance<AmorPhaseMaterial> ( *this );
}