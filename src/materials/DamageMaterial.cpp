/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic elastic-damage material
 *  This represents the material at a point in space.
 *  It is implemented in such a way that can be used for any
 *  discretisation methods, say finite elements, EFG and so on.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

/** @file DamageMaterial.cpp
 *  @brief Isotropic damage material model.
 * 
 *  Original Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *  
 *  Updates (when, what and who)
 *     - [21 February 2024] Integrated old damage laws,
 *     linear1, linear2, and new polynomial and multilinear.
 *     Removed Rankine type equivalent strain measure. (RB)
 */

#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/utilities.h>
#include <jem/io/NumberFormat.h>
#include <algorithm>

#include "util/Constants.h"
#include "util/Invariants.h"
#include "util/MathUtils.h"
#include "util/TupleUtils.h"
#include "DamageMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;
using jem::util::Properties;
using jive::Vector;


//=======================================================================
//   class DamageMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  DamageMaterial::SOFTENING_PROP        = "softening";
const char*  DamageMaterial::EQUISTRAIN_PROP       = "eqvStrain";
const char*  DamageMaterial::KAPPAI_PROP           = "kappaI";
const char*  DamageMaterial::KAPPAP_PROP           = "kappaP";
const char*  DamageMaterial::KAPPAC_PROP           = "kappaC";
const char*  DamageMaterial::TENSILE_PROP          = "ft";
const char*  DamageMaterial::FRACTURE_ENERGY_PROP  = "gf";
const char*  DamageMaterial::ALPHA_PROP            = "alpha";
const char*  DamageMaterial::BETA_PROP             = "beta";
const char*  DamageMaterial::ETA_PROP              = "eta";
const char*  DamageMaterial::B_PROP                = "b";
const char*  DamageMaterial::LENGTH_PROP           = "lengthScale";
const char*  DamageMaterial::REMOVE_DAMAGE_PROP    = "threshold";

const char*  DamageMaterial::CRACK_WIDTH_PROP      = "crackWidth";

const char*  DamageMaterial::MAZARS_EQUI_STRAIN    = "Mazars";
const char*  DamageMaterial::MISES_EQUI_STRAIN     = "vonMises";
const char*  DamageMaterial::RANKINE_EQUI_STRAIN   = "Rankine";
const char*  DamageMaterial::ENERGY_STRAIN         = "Energy";

const char*  DamageMaterial::PERFECT_SOFTENING     = "perfect";
const char*  DamageMaterial::LINEAR1_SOFTENING     = "linear1";
const char*  DamageMaterial::LINEAR2_SOFTENING     = "linear2";
const char*  DamageMaterial::EXPONENT1_SOFTENING   = "exponential1"; 
const char*  DamageMaterial::EXPONENT2_SOFTENING   = "exponential2"; 
const char*  DamageMaterial::EXPONENT3_SOFTENING   = "exponential2Reg";
const char*  DamageMaterial::POWER_SOFTENING       = "power";
const char*  DamageMaterial::HYPERBOLIC_SOFTENING  = "hyperbolic";
const char*  DamageMaterial::EXPOENERGY_SOFTENING  = "expoEnergy";
const char*  DamageMaterial::POLYNOMIAL_SOFTENING  = "polynomial";
const char*  DamageMaterial::MULTILINEAR_SOFTENING = "multilinear";

const double DamageMaterial::CRITICAL_DAMAGE      = 0.999;

vector<String> DamageMaterial::equiStrainDefs   (3);
vector<String> DamageMaterial::softeningLawDefs (11);


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


DamageMaterial::DamageMaterial

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  kappaI_     = 0.0;
  kappaP_     = 0.0;
  kappaC_     = 0.0;
  alpha_      = 0.0;
  beta_       = 0.0;
  eta_        = 0.0;
  b_          = 0.0;
  c_          = 0.0;
  threshold_  = 0.1;

  softening_  = Linear1;
  equiStrain_ = VonMises;

  elasticMat_ = newInstance<HookeMaterial> ( rank, globdat );

  elasticMod_ .resize ( STRAIN_COUNTS[rank], STRAIN_COUNTS[rank] );

  equiStrainDefs[0]   = MAZARS_EQUI_STRAIN;
  equiStrainDefs[1]   = MISES_EQUI_STRAIN;
  equiStrainDefs[2]   = ENERGY_STRAIN;

  softeningLawDefs[0] = PERFECT_SOFTENING;
  softeningLawDefs[1] = LINEAR1_SOFTENING;
  softeningLawDefs[2] = LINEAR2_SOFTENING;
  softeningLawDefs[3] = EXPONENT1_SOFTENING;
  softeningLawDefs[4] = EXPONENT2_SOFTENING;
  softeningLawDefs[5] = POWER_SOFTENING;
  softeningLawDefs[6] = HYPERBOLIC_SOFTENING;
  softeningLawDefs[7] = EXPONENT3_SOFTENING;
  softeningLawDefs[8] = EXPOENERGY_SOFTENING;   // German, IJNME 07
  softeningLawDefs[9] = POLYNOMIAL_SOFTENING;   // COMSOL
  softeningLawDefs[10] = MULTILINEAR_SOFTENING; // COMSOL
}


DamageMaterial::~DamageMaterial ()
{} 


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DamageMaterial::configure ( const Properties& props,
                                 const Properties& globdat )
{
  using jem::maxOf;
  using std::find;

  elasticMat_->configure ( props, globdat );

  bool ki = props.find ( kappaI_, KAPPAI_PROP, 0.0, maxOf( kappaI_) );

  if ( !ki )
  {
    props.get ( ft_, TENSILE_PROP, 0.0, maxOf ( ft_ ) );

    kappaI_ = ft_ / elasticMat_->getYoung ();
  }
  else
  {
    ft_ = kappaI_ * elasticMat_->getYoung ();
  }

  props.find ( gf_,     FRACTURE_ENERGY_PROP, 0.0, maxOf ( gf_) );
  props.find ( lambda_, CRACK_WIDTH_PROP,     0.0, 1000.0       );
  props.find ( c_,      LENGTH_PROP,          0.0, maxOf ( c_)  );

  // read type of softening law used

  String      softening;

  props.get ( softening, SOFTENING_PROP  );

  if ( find ( softeningLawDefs.begin (),
              softeningLawDefs.end   (),
              softening ) == softeningLawDefs.end () )
  {
    throw Error (
      JEM_FUNC,
      String("unexpected definition of softening law!!!\n") +
      String("Supported softening laws include: \n") +
        PERFECT_SOFTENING     + String(", ") 
      + LINEAR1_SOFTENING     + String(", ") 
      + LINEAR2_SOFTENING     + String(", ") 
      + EXPONENT1_SOFTENING   + String(", ") 
      + EXPONENT2_SOFTENING   + String(", ") 
      + POWER_SOFTENING       + String(", ") 
      + HYPERBOLIC_SOFTENING  + String(", ") 
      + POLYNOMIAL_SOFTENING  + String(", ") 
      + MULTILINEAR_SOFTENING
    );
  }

  // define correctly the parameters for each softening law

  if      ( softening == EXPONENT1_SOFTENING )
  {
    props.get ( alpha_,  ALPHA_PROP,  0.0, maxOf ( alpha_  ) );
    props.get ( beta_,   BETA_PROP,   0.0, maxOf ( beta_   ) );

    softening_ = Exponential3Params;
  }
  else if ( softening == EXPONENT2_SOFTENING )
  {
    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );

    softening_ = Exponential2Params;
  }
  else if ( softening == EXPONENT3_SOFTENING )
  {
    softening_ = Exponential2ParamsReg;
  }
  else if ( softening == POWER_SOFTENING )
  {
    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );
    props.get ( alpha_,  ALPHA_PROP,  0.0, maxOf ( alpha_  ) );
    props.get ( beta_,   BETA_PROP,   0.0, maxOf ( beta_   ) );

    softening_ = Power;
  }
  else if ( softening == HYPERBOLIC_SOFTENING )
  {
    props.get ( b_, B_PROP, 0.0, maxOf ( b_ ) );

    softening_ = Hyperbolic;
  }
  else if ( softening == LINEAR1_SOFTENING )
  {
    softening_ = Linear1;

    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );
  }
  else if ( softening == LINEAR2_SOFTENING )
  {
    softening_ = Linear2;

    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );
  }
  else if ( softening == EXPOENERGY_SOFTENING )
  {
    props.get ( kappaI_,  KAPPAI_PROP,  0.0, maxOf ( kappaI_  ) );
    props.get ( kappaC_,  KAPPAC_PROP,  0.0, maxOf ( kappaC_  ) );
    props.get ( s_,       "s",          0.0, maxOf ( s_       ) );

    softening_ = ExpoEnergy;
  }
  else if ( softening == POLYNOMIAL_SOFTENING )
  {
    softening_ = Polynomial;

    props.get ( kappaC_, KAPPAC_PROP, 0.0, maxOf ( kappaC_ ) );
  }
  else if ( softening == MULTILINEAR_SOFTENING )
  {
    softening_ = Multilinear;

    props.get ( kappaP_, KAPPAP_PROP, 0.0,     maxOf ( kappaP_ ) );
    props.get ( kappaC_, KAPPAC_PROP, kappaP_, maxOf ( kappaC_ ) );
  }

  // read type of equivalent strain definition

  String equiStrain;

  props.get ( equiStrain, EQUISTRAIN_PROP );

  if ( find ( equiStrainDefs.begin (), 
              equiStrainDefs.end   (), 
              equiStrain ) == equiStrainDefs.end () )
  {
    throw Error (
      JEM_FUNC,
      String ("unexpected definition of equivalent strain!!!\n") +
      String ("Available ones are: ") 
      + MAZARS_EQUI_STRAIN + String(", ") 
      + MISES_EQUI_STRAIN  + String(", ") 
      + ENERGY_STRAIN
    );
  }

  if      ( equiStrain == MISES_EQUI_STRAIN )
  {
    props.get ( eta_, ETA_PROP, 0.0, maxOf ( eta_ ) );

    equiStrain_ = VonMises;
  }
  else if ( equiStrain == MAZARS_EQUI_STRAIN )
  {
    equiStrain_ = Mazars;
  }
  else if ( equiStrain == ENERGY_STRAIN )
  {
    equiStrain_ = StrainEnergy;
  }

  vonMisesEqv_.init ( eta_, elasticMat_->getPoisson ( ) );

  elasticMod_ = elasticMat_->getStiffMat();
    
  props.find ( threshold_, REMOVE_DAMAGE_PROP, 0.0, 1.0 );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DamageMaterial::getConfig ( const Properties& conf,
                                 const Properties& globdat ) const
{
  elasticMat_->getConfig ( conf , globdat);

  conf.set ( SOFTENING_PROP,  softening_  );
  conf.set ( EQUISTRAIN_PROP, equiStrain_ );

  conf.set ( KAPPAI_PROP, kappaI_ );
  conf.set ( KAPPAC_PROP, kappaC_ );
  conf.set ( ALPHA_PROP,  alpha_  );
  conf.set ( BETA_PROP,   beta_   );
  conf.set ( ETA_PROP,    eta_    );
  conf.set ( REMOVE_DAMAGE_PROP, threshold_ );

  if ( c_ != 0 ) conf.set ( LENGTH_PROP, c_ );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void DamageMaterial::updateConfig ()

{
  elasticMat_-> updateConfig();
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  DamageMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )
{
  // Get the nonlocal equivalent strain:

  const int n          = stress.size ();
  const int m          = strain.size ();

  double    equiStrain;

  if ( n == m ) // local damage models
  {
    equiStrain = getEquiStrain ( strain );
  }
  else
  {
    equiStrain = strain[n];
  }

  // Get the history variable of the previous converged load step

  double    kappa0     = preHist_.eqveps[ipoint];

  // Compute the loading function f

  double    f          = equiStrain - kappa0;

  double kappa;
  int    load;

  // f == 0 (exactly equals to zero), at the beginning of every 
  // loading step, this condition holds for loading integration point

  if ( f == 0.0 )
  {
    kappa = kappa0;
    load  = preHist_.loading[ipoint];
  }  
  else 
  {
    // f > 0: loading, update history variable

    if ( f > 0.0 )
    {  
      kappa = equiStrain;  
      load  = ( kappa < kappaI_ ) ? 0 : 1;
    }
 
    // f < 0: unloading, keep the same history variable

    else 
    {
      kappa = kappa0;
      load  = 0;
    }
  }

  // compute the damage variable

  double omega  = damageEvolution_ ( kappa );

  // secant stiffness matrix only
  // the second term will be computed by non-local integral  model

  stiff         = ( 1.0 - omega ) * elasticMod_;

  // Compute stress vector

  MatmulChain<double,1>     mc1;

  stress = mc1.matmul ( stiff, strain[slice(BEGIN,n)] );

  // update history variables

  newHist_.eqveps[ipoint]  = kappa;
  newHist_.loading[ipoint] = load;

}

// ----------------------------------------------------------------
//  update (overloaded version)
// ----------------------------------------------------------------

void  DamageMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint,
      double                he )
{
  MatmulChain<double,1>     mc1;
  MatmulChain<double,2>     mc2;

  double    equiStrain = getEquiStrain ( strain ); // Get the local equivalent strain:
  double    kappa0     = preHist_.eqveps[ipoint];  // Get the history variable of the previous converged load step
  double    f          = equiStrain - kappa0;      // Compute the loading function f

  double kappa;
  int    load;

  // f == 0 (exactly equals to zero), at the beginning of every 
  // loading step, this condition holds for loading integration point

  if ( f == 0.0 )
  {
    kappa = kappa0;
    load  = preHist_.loading[ipoint];
  }  
  else 
  {
    if ( f > 0.0 )
    {
      kappa = equiStrain;  
      load  = ( kappa < kappaI_ ) ? 0 : 1;
    }
 
    // f < 0: unloading, keep the same history variable

    else 
    {
      kappa = kappa0;
      load  = 0;
    }
  }

  // compute the damage variable

  double omega  = damageEvolution_ ( kappa );

  // secant stiffness matrix only
  // the second term will be computed by non-local integral  model

  stiff         = ( 1.0 - omega ) * elasticMod_;
        
  if ( load )
  {
     double  dOmega      = getdOmegadKappa ( kappa   );
     Vector  dEpsBardEps = getdEpsBardEps  ( strain  );	
  
     stiff -= dOmega * mc2.matmul(elasticMod_,matmul(strain,dEpsBardEps));	
  }

  // Compute stress vector

  stress = ( 1.0 - omega ) * mc1.matmul ( elasticMod_, strain );

  // update history variables

  newHist_.eqveps[ipoint]  = kappa;
  newHist_.loading[ipoint] = load;

}

// ----------------------------------------------------------------
//  update (overloaded version)
// ----------------------------------------------------------------

void  DamageMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain0,
      const Vector&         dstrain,
      int                   ipoint,
      double                 he )
{
    Vector tstrain (STRAIN_COUNTS[rank_]);
    tstrain = strain0 + dstrain;
    update ( stress, stiff, tstrain, ipoint, he );
}


//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void DamageMaterial::allocPoints ( int count )
{
  preHist_.eqveps. resize ( count );
  preHist_.loading.resize ( count );

  newHist_.eqveps. resize ( count );
  newHist_.loading.resize ( count );

  preHist_.eqveps  = 0.0;
  preHist_.loading = 0;

  newHist_.eqveps  = 0.0;
  newHist_.loading = 0;
}


// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  DamageMaterial::commit()
{
  newHist_.eqveps. swap ( preHist_.eqveps  );
  newHist_.loading.swap ( preHist_.loading );
}


//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> DamageMaterial::clone ( ) const

{
  return newInstance<DamageMaterial> ( *this );
}


// ---------------------------------------------------------------
//  isFullyDamaged ( returns true for fully damaged ipoint )
// ---------------------------------------------------------------

bool DamageMaterial::isFullyDamaged ( int ipoint ) const
{
  double kappa = newHist_.eqveps[ ipoint ];
  double omega = damageEvolution_ ( kappa );

  return ( 1.0 - omega < threshold_ ) ? true : false ;
}


// -------------------------------------------------------------------
//  giveElasticMatrix
// -------------------------------------------------------------------

Matrix  DamageMaterial::giveElasticMatrix ( ) const
{
  return elasticMod_;
}


//-----------------------------------------------------------------------
//   getEquiStrain
//-----------------------------------------------------------------------

double DamageMaterial::getEquiStrain

  ( const Vector& strain ) 

{
  double eqvStr = 0.0;

  if      ( equiStrain_ == VonMises  )
  {
    Invariants strInvs ( strain, /*isStrain =*/ true );

    const double I1 = strInvs.getI1 ();
    const double J2 = strInvs.getJ2 ();

    eqvStr     = vonMisesEqv_ ( I1, J2 );
  } 
  else if ( equiStrain_ == Mazars  )
  {
    if ( strain.size() == 1 )
    {
      return  ::fabs ( strain[0] );
    }

    Invariants strInvs ( strain, /*isStrain =*/ true );

    Vec3 princ = strInvs.getPrincipalValues ();

    eqvStr  = 0.0;

    for ( int i = 0; i < 3; i++ )
    {
      double x =  princ[i] ;

      if ( x > 0 )
      {
      eqvStr += x * x ;
      }
    }

    eqvStr = ::sqrt ( eqvStr );
  }

  else if ( equiStrain_ == StrainEnergy  )
  {
    MatmulChain<double,1>     mc1;

    eqvStr        = 0.5 * dot ( strain, mc1.matmul ( elasticMod_, strain ) );
  }
 
  return eqvStr;
}

// --------------------------------------------------------------------
//  getdEpsBardEps (derivatives of equivalent strain w.r.t strain vector)
// --------------------------------------------------------------------

Vector DamageMaterial::getdEpsBardEps

  ( const Vector& strain ) 

{
  const int strCount = strain.size ( );

  Vector  ret ( strCount );

  if      ( equiStrain_ == VonMises  )
  {
    Invariants strInvs ( strain, /*isStrain =*/ true );

    const double I1   = strInvs.getI1 ();
    const double J2   = strInvs.getJ2 ();

    const Vec6 GradI1 = strInvs.getGradI1 ();
    const Vec6 GradJ2 = strInvs.getGradJ2 ();

    Vector dI1 (6);
    Vector dJ2 (6);

    t6_to_v6( dI1, GradI1 );
    t6_to_v6( dJ2, GradJ2 );

    ret  = vonMisesEqv_ ( dI1, dJ2, I1, J2 );
  }

  else if ( equiStrain_ == Mazars  )
  {
    if ( strCount == 1 ) 
    {
      ret = ( strain[0] > 0.0 )? 1.0 : 0.0; 
    }
    else
    {
      ret = getdEquidEpsMazars_ ( strain );
    }
  } 

  else if ( equiStrain_ == StrainEnergy  )
  {
    MatmulChain<double,1>     mc1;

    ret = mc1.matmul ( elasticMod_, strain );;
  }

  return ret;
}

//-----------------------------------------------------------------------
//   getdOmegadKappa
//-----------------------------------------------------------------------

double DamageMaterial:: getdOmegadKappa

    ( double k )  const

{
  double dOmega = 0.0;

  if           ( softening_ == Perfect  )
  {
    dOmega = getDerivPerfectSoftening( kappaI_, k );
  } 
  else if      ( softening_ ==  Linear1 )
  {
    dOmega = getDerivLinear1Softening( kappaI_, kappaC_ , k );
  }
  else if      ( softening_ ==  Linear2 )
  {
    dOmega = getDerivLinear2Softening( kappaI_, kappaC_ , k );
  } 
  else if      ( softening_ == Exponential3Params )
  {
    dOmega = getDerivExponentialSoftening1 ( kappaI_, alpha_, beta_, k );
  }
  else if      ( softening_ == Exponential2Params )
  {
    dOmega = getDerivExponentialSoftening2 ( kappaI_, kappaC_, k );
  }
  else if      ( softening_ == Power )
  {
    dOmega = getDerivPowerSoftening ( kappaI_, kappaC_, alpha_, beta_, k );
  }
  else if      ( softening_ == Hyperbolic )
  {
    dOmega = getDerivHyperbolicSoftening ( kappaI_, b_,k );
  }
  else if      ( softening_ == ExpoEnergy )
  {
    dOmega = getDerivExpoEnergySoftening ( kappaI_, kappaC_, s_, k );
  }
  else if      ( softening_ ==  Polynomial )
  {
    dOmega = getDerivPolynomialSoftening( kappaI_, kappaC_, k );
  }
  else if      ( softening_ ==  Multilinear )
  {
    dOmega = getDerivMultilinearSoftening( kappaI_, kappaP_, kappaC_, k );
  }

  return dOmega;
}

//-----------------------------------------------------------------------
//   getdOmegadKappa (for regularized local damage)
//-----------------------------------------------------------------------

double DamageMaterial:: getdOmegadKappa

    ( double k, double he )  const

{
  double dOmega = 0.0;

  if       ( softening_ == Exponential3Params )
  {
    dOmega = getDerivExponentialSoftening1 ( kappaI_, alpha_, gf_, ft_, k, he );
  }
  else if  ( softening_ == Exponential2Params )
  {
    dOmega = getDerivExponentialSoftening2 
                ( kappaI_, kappaC_, k, lambda_, he );
  }
  else if  ( softening_ == Exponential2ParamsReg )
  {
    dOmega = getDerivExponentialSoftening3
                ( kappaI_, k, gf_, ft_, he );
  }

  return dOmega;
}


//-----------------------------------------------------------------------
//   damageEvolution_
//-----------------------------------------------------------------------

// Remark: The following Fortran-like implementation is ugly
// but efficient. There is an alternative using policy-based
// class by writing orthogonal classes for every softening
// law. It is, however, slower than the current implementation.

double DamageMaterial::damageEvolution_

  ( double k ) const

{
  double omega = 0.0;

  if      ( softening_ == Perfect )
  {
    omega = perfectSoftening( kappaI_, k ); 
  }
  else if ( softening_ == Linear1 )
  {
    omega = linear1Softening( kappaI_, kappaC_ , k );
  }
  else if ( softening_ == Linear2 )
  {
    omega = linear2Softening( kappaI_, kappaC_ , k );
  } 
  else if ( softening_ == Exponential3Params )
  {
    omega = exponentialSoftening1 ( kappaI_, alpha_, beta_, k );
  }
  else if ( softening_ == Exponential2Params )
  {
    omega = exponentialSoftening2 ( kappaI_, kappaC_ , k );
  }
  else if ( softening_ == Power )
  {
    omega = powerSoftening ( kappaI_, kappaC_, alpha_, beta_, k );
  }
  else if ( softening_ == Hyperbolic )
  {
    omega = hyperbolicSoftening ( kappaI_, b_, k );
  }
  else if ( softening_ == ExpoEnergy )
  {
    omega = expoEnergySoftening ( kappaI_, kappaC_, s_, k );
  }
  else if ( softening_ == Polynomial )
  {
    omega = polynomialSoftening( kappaI_, kappaC_, k );
  }
  else if ( softening_ == Multilinear )
  {
    omega = multilinearSoftening( kappaI_, kappaP_, kappaC_, k );
  }

  return omega;
}

//-----------------------------------------------------------------------
//   damageEvolution_
//-----------------------------------------------------------------------

double DamageMaterial::damageEvolution_

  ( double k, double he ) const

{
  double omega = 0.0;

  if      ( softening_ == Exponential3Params )
  {
    omega = exponentialSoftening1 ( kappaI_, alpha_, gf_, ft_ , k, he );
  }
  else if ( softening_ == Exponential2Params )
  {
    omega = exponentialSoftening2 
             ( kappaI_, kappaC_ , k, lambda_, he );
  }
  else if ( softening_ == Exponential2ParamsReg )
  {
    omega = exponentialSoftening3 ( kappaI_, k, gf_, ft_, he );
  }

  return omega;
}


// ---------------------------------------------------------------
// compute derivatives of Mazars equivalent strain w.r.t strain
// ---------------------------------------------------------------

Vector  DamageMaterial:: getdEquidEpsMazars_

  ( const Vector&   strain )

{
  const int s = strain.size();

  Vector   ret ( s );

  ret   = 0.0 ;

  // 3D problem

  if ( s == 6 )
  {
    Invariants inv ( strain, /*isStrain =*/ true );

    const Vec3 princ  = inv.getPrincipalValues ();

    const double I1  = inv.getI1 ();
    const double I2  = inv.getI2 ();

    const Vec6 gradI1 = inv.getGradI1 ();
    const Vec6 gradI2 = inv.getGradI2 ();
    const Vec6 gradI3 = inv.getGradI3 ();

    Vector dI1 (6);
    Vector dI2 (6);
    Vector dI3 (6);

    t6_to_v6( dI1, gradI1 );
    t6_to_v6( dI2, gradI2 );
    t6_to_v6( dI3, gradI3 );

    Vector    dEpsIdEps (6);

    double epsBar = 0.0;

    for ( int i = 0; i < 3; i++ )
    {
      double x =  princ[i] ;

      if ( x > 0 )
      {
	    epsBar += x * x ;
      }
    }

    epsBar = ::sqrt ( epsBar );

    // check for case of numerically zero strain vector

    if ( epsBar < EPS )
    {
      return ret; 
    }

    double vi, temp;
 
    for ( int i = 0 ; i < 3; i++ )
    {
      vi  = princ[i];

      if ( vi > 0 )
      {
	    temp      = 1.0 / ( 3.0 * vi * vi - 2.0 * I1 * vi + I2);
	    dEpsIdEps = temp * ( vi * vi * dI1 - vi * dI2 + dI3 );
	   ret       += vi * dEpsIdEps ;
      }
    }

    ret          *= (1.0 / epsBar) ;
  }
  else
  {
    // 2D: strain={eps_xx, eps_yy, eps_zz, eps_xy}  

    double exx =       strain[0];
    double eyy =       strain[1];
    double exy = 0.5 * strain[3];

    double exxMeyy = exx - eyy;
    double exxPeyy = exx + eyy;

    double poi  = elasticMat_->getPoisson();
    double poiD = poi / ( poi - 1.0 );

    // principal strains are roots of a quadratic equation with det = d

    double d   = ::sqrt ( exxMeyy * exxMeyy + 4.0 * exy * exy );

    double prinstr0 = 0.5 * ( exxPeyy + d );
    double prinstr1 = exxPeyy - prinstr0;
    double prinstr2 = elasticMat_-> getState() == PlaneStress ? poiD * exxPeyy : 0.0;

    prinstr0 = prinstr0 < 0. ? 0. : prinstr0;
    prinstr1 = prinstr1 < 0. ? 0. : prinstr1;
    prinstr2 = prinstr2 < 0. ? 0. : prinstr2;

    double den = ::sqrt ( prinstr0 * prinstr0 + prinstr1 * prinstr1 + prinstr2 * prinstr2 ); 

    // check for case of numericallly zero strain vector

    if ( den < 1e-08 )
    {
      return ret; 
    }

    // derivatives of equivalent strain w.r.t principal strains

    double denInv = 1. / den;

    double dedprin0 = prinstr0 * denInv;
    double dedprin1 = prinstr1 * denInv;
    double dedprin2 = prinstr2 * denInv;

    // derivatives of principal strains w.r.t strain tensor
    // avoid division by zero

    if ( d == 0.0 )
    {
      d = 1.0;
    }

    double fac     =  0.5 / d; 

    double de1dexx = 0.5 + fac * exxMeyy;
    double de1deyy = 1.0 - de1dexx;
    double de1dexy = 2.0 * fac * exy ; // attention here, 2 = 4 (from derivation) * 0.5

    double de2dexx =  de1deyy;
    double de2deyy =  de1dexx;
    double de2dexy = -de1dexy;

    double de3dexx = poiD;
    double de3deyy = de3dexx;

    // finally, derivatives of equivalent strain w.r.t strain tensor
    // by chain rule

    ret[0] = dedprin0 * de1dexx +  dedprin1 * de2dexx +  dedprin2 * de3dexx;
    ret[1] = dedprin0 * de1deyy +  dedprin1 * de2deyy +  dedprin2 * de3deyy;
    ret[3] = dedprin0 * de1dexy +  dedprin1 * de2dexy;
  }

  return ret;
}


// ======================================================================
//   Implementation of related functions
// ======================================================================


//-----------------------------------------------------------------------
//   perfectSoftening
//-----------------------------------------------------------------------

double                      perfectSoftening

  ( double Ki, double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = 1.0 - Ki / K;

  return omega;
}

//-----------------------------------------------------------------------
//   getDerivPerfectSoftening
//-----------------------------------------------------------------------

double                      getDerivPerfectSoftening

    ( double Ki, double K )

{
  if ( K < Ki )
  {
    return 0.0;
  }

  return  Ki / ( K * K );
}


//-----------------------------------------------------------------------
//   linearSoftening (type I)
//-----------------------------------------------------------------------

double                      linear1Softening

  ( double Ki,
    double Kc,
    double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = K > Kc ? 1.0 : Kc * (K - Ki) / (K * (Kc - Ki));
 
  if ( omega == 1.0 )
  {
    omega = DamageMaterial::CRITICAL_DAMAGE;
  }

  return omega;
}

//-----------------------------------------------------------------------
//   getDerivLinearSoftening (type I)
//-----------------------------------------------------------------------

double                      getDerivLinear1Softening

    ( double Ki,
      double Kc,
      double K )

{
  return  (K < Kc) ?  ( Ki * Kc ) / ( (Kc - Ki) * K * K ) : 0.0;
}


//-----------------------------------------------------------------------
//   linear2Softening
//-----------------------------------------------------------------------

double                      linear2Softening

  ( const double Ki,
    const double Kc,
    const double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }
  else
  {
    return K > Kc ? DamageMaterial::CRITICAL_DAMAGE 
                  : (K - Ki) / (Kc - Ki);
  }
}


//-----------------------------------------------------------------------
//   getDerivLinear2Softening
//-----------------------------------------------------------------------

double                      getDerivLinear2Softening

    ( const double Ki,
      const double Kc,
      const double K )

{
  return  (K < Kc) ?  1.0 / (Kc - Ki) : 0.0;
}


//-----------------------------------------------------------------------
//   exponentialSoftening1
//-----------------------------------------------------------------------

double                      exponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    double K )
{
  if ( K < Ki ) return 0.0;

  double omega = 1.0 - (Ki / K) * ( 1.0 - alpha + alpha * exp (-beta * (K - Ki) ) );

  //return isTiny ( 1.0 - omega ) ? DamageMaterial::CRITICAL_DAMAGE : omega;

  return omega;
}

// overloaded version

double                      exponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K,
    double he )
{
  if ( K < Ki ) return 0.0;

  double beta  = alpha / ( -0.5 * Ki + gf / ft /he );

  double omega = 1.0 - (Ki / K) * ( 1.0 - alpha + alpha * exp (-beta * (K - Ki) ) );

  return omega;
}

//-----------------------------------------------------------------------
//   getDerivExponentialSoftening1
//-----------------------------------------------------------------------

double                      getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    double K )
{
  const double h = K - Ki;	
  const double a = Ki / (K * K);
  const double b = a * K;
  const double c = alpha * exp(-beta * h );

  return a * (1.0 - alpha + c) + beta * b * c;
}

// overloaded version

double                      getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K,
    double he )
{
  double beta    = alpha / ( -0.5 * Ki + gf / ft /he );

  const double h = K - Ki;	
  const double a = Ki / (K * K);
  const double b = a * K;
  const double c = alpha * exp(-beta * h );

  return a * (1.0 - alpha + c) + beta * b * c;
}

//-----------------------------------------------------------------------
//    exponentialSoftening2
//-----------------------------------------------------------------------

double                       exponentialSoftening2

  ( double Ki,
    double Kc,
    double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = 1.0 - (Ki / K) * exp ( (Ki - K) / (Kc - Ki) );

  return omega;
}

// overloaded version

double                       exponentialSoftening2

  ( double Ki,
    double Kc,
    double K, double lambda,
    double he )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double KC    = lambda / he * ( Kc - 0.5 * Ki ) + 0.5 * Ki;

  double omega = 1.0 - (Ki / K) * exp ( (Ki - K) / (KC - Ki) );

  return omega;  
}

//-----------------------------------------------------------------------
//    getDerivExponentialSoftening2
//-----------------------------------------------------------------------

double                       getDerivExponentialSoftening2

    ( double Ki,
      double Kc,
      double K )

{ 
  double a = Ki - K;
  double b = Kc - Ki;
  double c = exp ( a / b );
  double d = Ki / (K * K);

  return d * c + d * K * c / b;
}

// overloaded version

double                       getDerivExponentialSoftening2

    ( double Ki,
      double Kc,
      double K, double lambda,
      double he )

{
  double KC= lambda / he * ( Kc - 0.5 * Ki ) + 0.5 * Ki;

  double a = Ki - K;
  double b = KC - Ki;
  double c = exp ( a / b );
  double d = Ki / (K * K);

  return d * c + d * K * c / b;
}

// Jirasek regularised exponential softening law

double                  exponentialSoftening3

  ( double Ki, double K,
    double gf, double ft, 
    double he )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double Kc    = 0.5 * Ki + gf / ft / he;
  double omega = 1.0 - (Ki / K) * exp ( (Ki - K) / (Kc - Ki) );

  return  omega;
}

double                  getDerivExponentialSoftening3

   ( double Ki, double K,
     double gf, double ft,
     double he )
{
  double Kc = 0.5 * Ki + gf / ft / he;

  double a  = Ki - K;
  double b  = Kc - Ki;
  double c  = exp ( a / b );
  double d  = Ki / (K * K);

  return d * c + d * K * c / b;
}

//-----------------------------------------------------------------------
//    hyperbolicSoftening
//-----------------------------------------------------------------------

double                       hyperbolicSoftening

  ( double Ki,
    double b,
    double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }

  double omega = 1.0 - 1.0 / ( 1.0 + b * ( K - Ki ) );

  return isTiny ( 1.0 - omega ) ? DamageMaterial::CRITICAL_DAMAGE : omega;  
}

//-----------------------------------------------------------------------
//    getDerivExponentialSoftening2
//-----------------------------------------------------------------------

double                       getDerivHyperbolicSoftening

    ( double Ki,
      double b,
      double K )

{

  return b / ( 1.0 + b * ( K - Ki ) ) / ( 1.0 + b * ( K - Ki ) ) ;
 
}


//-----------------------------------------------------------------------
//  powerSoftening
//-----------------------------------------------------------------------

// particularly implemented for the composite compact tension test

double                      powerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K )
{
  double omega;

  if ( K < Ki )
  {
    omega = 0.0;
  }
  else
  {
    double Kcc    =  Kc;
    double omegac = 1.0 - ::pow ( Ki / Kcc, beta ) * 
                          ::pow ( (Kc - Kcc) / (Kc - Ki) , alpha );

    omega  = (K > Kcc) ? omegac : 1.0 -::pow ( Ki / K, beta ) * 
                                       ::pow ( (Kc - K) / (Kc - Ki) , alpha );
  }

  return omega; 

  /*

  // the damage evolution in Geers - Enhanced solution control

  double omega;

  if ( K < Ki )
  {
    omega = 0.0;
  }
  else
  {
    const double gamma = 0.01;
    const double eta   = 1.0 - gamma;
    const int    n     = -5;

    double  muy  = ( gamma - n * eta ) / ( n * eta * ::pow ( Ki, n ) );

    omega        = 1.0 - gamma * Ki / K - eta * ::pow ( K / Ki, n) *
                    ::exp ( muy * ( ::pow ( K, n ) - ::pow( Ki, n ) ) );  
  }
  
  return omega; 
  */
  
}

//-----------------------------------------------------------------------
//  getDerivPowerSoftening
//-----------------------------------------------------------------------

double                      getDerivPowerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K )
{
  double ret;
  double Kcc       =  Kc;

  if ( K > Kcc )
  {
    ret =  0.0;
  }
  else
  {
    const double a = Ki / K ;
    const double b = 1.0 / ( Kc - Ki );
    const double c = ( Kc - K ) * b;
    const double d = ::pow ( a, beta  - 1.0 );
    const double e = ::pow ( c, alpha - 1.0 );

    ret  = d * e * ( (a / K) * beta * c  +  a * alpha * b );
  }

  return ret ; 

/*
  const double gamma = 0.01;
  const double eta   = 1.0 - gamma;
  const int    n     = -5;

  double  muy  = ( gamma - n * eta ) / ( n * eta * ::pow ( Ki, n ) );

  double fac1  =  exp ( muy * ( ::pow ( K, n )  - ::pow( Ki, n ) ) );

  double ret   =  gamma * Ki / K / K - eta * n * ::pow ( K / Ki, n-1) / Ki * fac1
                 - eta * pow ( K / Ki, n) * fac1 * muy * n * ::pow ( K, n-1);
 
  return ret;
*/
}

//-----------------------------------------------------------------------
//   expoEnergySoftening
//-----------------------------------------------------------------------


double                     expoEnergySoftening

 ( double kappaI, double kappaC,
   double s, double kappa )
{
  if ( kappa < kappaI ) return 0.;

  double tem = pow ( ( kappa - kappaI ) / kappaC, s );

  return ( 1. - exp ( - tem ) );
}

//-----------------------------------------------------------------------
//   getDerivExpoEnergySoftening
//-----------------------------------------------------------------------


double                      getDerivExpoEnergySoftening

  ( double kappaI, double kappaC,
    double s, double kappa )
{
  if ( kappa < kappaI ) return 0.;

  double tem1 = pow ( ( kappa - kappaI ) / kappaC, s   );
  double tem2 = pow ( ( kappa - kappaI ) / kappaC, s-1 );

  return  exp ( - tem1 ) * tem2 * s / kappaC;
}


//-----------------------------------------------------------------------
//   polynomialSoftening
//-----------------------------------------------------------------------

double                      polynomialSoftening

  ( const double Ki,
    const double Kc,
    const double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }
  else
  {
    const double a = (K - Ki) / (Kc - Ki);

    return 1.0 - (Ki / K) * ( 1. + ( ::pow(a, 2.) * (2. * a - 3.) ) );
  }
}


//-----------------------------------------------------------------------
//   getDerivPolynomialSoftening
//-----------------------------------------------------------------------

double                      getDerivPolynomialSoftening

    ( const double Ki,
      const double Kc,
      const double K )

{
  if ( K < Ki )
  {
    return 0.0;
  }
  else
  {
    const double a = (Ki - K);
    const double b = (Kc - Ki);
    const double c = a/b;
    const double d = 1. + ( ::pow(a, 2.) * (2. * a - 3.) );

    return - Ki/(K*K) * (  2. * K * ( c/b * ( 2. * c -3. ) + ::pow(c,2.) / b ) 
                         - d );
  }
}


//-----------------------------------------------------------------------
//   multilinearSoftening
//-----------------------------------------------------------------------

double                      multilinearSoftening

  ( const double Ki,
    const double Kp,
    const double Kc,
    const double K )
{
  if ( K < Ki )
  {
    return 0.0;
  }
  else 
  {
    if ( K > Kp )
    {
      return perfectSoftening (Ki, K);
    }
    else
    {
      return linear2Softening (Kp, Kc, K);
    }
  }
}


//-----------------------------------------------------------------------
//   getDerivMultilinearSoftening
//-----------------------------------------------------------------------

double                      getDerivMultilinearSoftening

    ( const double Ki,
      const double Kp,
      const double Kc,
      const double K )

{
  if ( K < Ki )
  {
    return 0.0;
  }
  else 
  {
    if ( K > Kp )
    {
      return getDerivPerfectSoftening (Ki, K);
    }
    else
    {
      return getDerivLinear2Softening (Kp, Kc, K);
    }
  }
}