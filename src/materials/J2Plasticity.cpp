
/** @file J2Plasticity.h
 *  @brief Implements a mixed isotropic-kinemetic linear hardening J2
 *  plasticity model.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 12 June 2023
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2023] 
 */

#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/util/Flex.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/utilities.h>

#include "util/BasicUtils.h"
#include "J2Plasticity.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  J2Plasticity::YOUNG_PROP          = "young";
const char*  J2Plasticity::POISSON_PROP        = "poisson";
const char*  J2Plasticity::YIELD_STRESS_PROP   = "sigY";
const char*  J2Plasticity::H_MODULUS_PROP      = "H";
const char*  J2Plasticity::R_PROP              = "r";
const char*  J2Plasticity::STATE_PROP          = "state";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


J2Plasticity::J2Plasticity 

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  // Compute elastic stiffness matrix

  elasticMat_ = newInstance<HookeMaterial> ( rank, globdat );
  
  elasticStiffMat_.resize( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  elasticStiffMat_ = 0.0;

  // Compute volumetric and deviatoric projection matrices

  P_.resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  P_ = 0.0;

  Q_.resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  Q_ = 0.0;

  for ( int i = 0; i<3; i++)
  {
    for ( int j = 0; j<3; j++)
    {
      i == j ? P_(i,j) = 2.0/3.0: P_(i,j) = -1.0/3.0;
      i == j ? Q_(i,j) = 2.0/3.0: Q_(i,j) = -1.0/3.0;
    }
  }

  if ( rank_ == 2 )
  {
    P_(3,3) = 2.0;
    Q_(3,3) = 0.5;
  }
  else if ( rank_ == 3 )
  {
    P_(3,3) = 2.0;
    P_(4,4) = 2.0;
    P_(5,5) = 2.0;

    Q_(3,3) = 0.5;
    Q_(4,4) = 0.5;
    Q_(5,5) = 0.5;
  }

}


J2Plasticity::~J2Plasticity ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void J2Plasticity::configure ( const Properties& props,
                               const Properties& globdat )
{
  props.find ( sigY_, YIELD_STRESS_PROP );
  props.find ( H_,    H_MODULUS_PROP    );
  props.find ( r_,    R_PROP, 0.0, 1.0  );

  elasticMat_->configure ( props, globdat );
  
  elasticStiffMat_ = elasticMat_->getStiffMat();

  G_ = young_ / 2. / ( 1. + poisson_ );
  K_ = young_ / 3. / ( 1. - 2. * poisson_ );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void J2Plasticity::getConfig ( const Properties& conf, 
                                const Properties& globdat ) const
{
  conf.set ( YIELD_STRESS_PROP, sigY_ );
  conf.set ( H_MODULUS_PROP   , H_    );
  conf.set ( R_PROP           , r_    );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void J2Plasticity::updateConfig ()

{
  // If young_ and poisson_ are updated with setters, they must be updated
  // in elasticMat_, followed by a call to compute the updated material
  // stiffness matrix.

  elasticStiffMat_ = elasticMat_->getStiffMat();

  G_ = young_ / 2. / ( 1. + poisson_ );
  K_ = young_ / 3. / ( 1. - 2. * poisson_ );
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void J2Plasticity::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{
  // NOTE: 'strain' refers to the strain increment over the last converged
  //       step, and not the total strain in the current step

  // Initialize a temporary stress variable (to be used for trial stress and
  // all subsequent manipulations.)

  Vector sigX ( STRAIN_COUNTS[rank_] );   sigX = 0.0;

  // Get old step values ( stress, back stress, kappa )

  Vector sig0 ( STRAIN_COUNTS[rank_] );   sig0 = preHist_[ipoint].sig;
  Vector alf0 ( STRAIN_COUNTS[rank_] );   alf0 = preHist_[ipoint].alf;
  double kap0;                            kap0 = preHist_[ipoint].kap;

  // Compute the reduced trial stress ( elastic step )

  sigX = sig0 + matmul ( elasticStiffMat_, strain ) - alf0;

  // Compute the equivalent stress and the yield function

  const double sigXEq = sqrt( std::abs( 3./2.*dot( sigX,matmul( P_,sigX ))));
  const double phiTr  = sigXEq - ( sigY_ + kap0 );

  // Check whether current stress is outside the yield surface

  if ( abs ( phiTr ) < 0.0 )
  {
    // Elastic step

    stress = sigX;
    stiff  = elasticStiffMat_;

    // Update internal variables

    newHist_[ipoint].sig     = sigX;
    newHist_[ipoint].epsp    = 0.0;
    newHist_[ipoint].kap     = kap0;
    newHist_[ipoint].alf     = alf0;
    newHist_[ipoint].loading = false;
  }
  else
  {
    // Plastic step

    // Initialize variables

    Vector f   ( STRAIN_COUNTS[rank_] );    f   = 0.0;
    Vector alf ( STRAIN_COUNTS[rank_] );    alf = 0.0;
    double kap;                             kap = 0.0;

    // Compute the integrated plastic multiplier

    const double mu = phiTr / ( 3. * G_ + H_ );

    // Compute f, sig, kap, alf

    f       = ( ( 3./2. ) / sigXEq ) * matmul ( P_, sigX );
    stress  = sigX - mu * 2. * G_ * f;
    kap     = kap0 + r_ * H_ * mu;
    alf     = alf0 + 2./3. * ( 1. - r_ ) * H_ * mu * f;

    // Compute algorithmic material tangent stiffness

    stiff = elasticStiffMat_ - 2. * G_ * ( 2. * G_ / ( 3. * G_ + H_ ) * ( kap0 + sigY_ ) / sigXEq
                                           * matmul ( f, f ) + ( mu * 3. * G_ / sigXEq ) * Q_ );

    // Update internal variables

    newHist_[ipoint].sig     = stress;
    newHist_[ipoint].epsp    = mu * f;
    newHist_[ipoint].kap     = kap;
    newHist_[ipoint].alf     = alf;
    newHist_[ipoint].loading = true;
  }
}


//-----------------------------------------------------------------------
//  allocPoints
//-----------------------------------------------------------------------

void J2Plasticity::allocPoints
    
  ( const idx_t          npoints )

{
  // Initialize history class

  Hist_ hist = Hist_();
  
  hist.sig .resize ( STRAIN_COUNTS[rank_] );
  hist.epsp.resize ( STRAIN_COUNTS[rank_] );
  hist.alf .resize ( STRAIN_COUNTS[rank_] );

  hist.sig      = 0.0;
  hist.epsp     = 0.0;
  hist.alf      = 0.0;
  hist.kap      = 0.0;
  hist.loading  = false;

  // Allocate hist for all integration points

  preHist_.reserve ( npoints );
  newHist_.reserve ( npoints );

  for ( idx_t ip = 0; ip < npoints; ++ip )
  {
    preHist_.pushBack ( hist );
    newHist_.pushBack ( hist );
  }

  latestHist_ = &preHist_;
}


//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> J2Plasticity::clone ( ) const

{
  return newInstance<J2Plasticity> ( *this );
}


//-----------------------------------------------------------------------
//  commit
//-----------------------------------------------------------------------


void J2Plasticity::commit ( )

{
  ( newHist_ ).swap ( preHist_ );

  latestHist_ =    &( preHist_ );
  
}


//-----------------------------------------------------------------------
//  getHistory
//-----------------------------------------------------------------------


void J2Plasticity::getHistory

  ( Vector&              hvals,
    const idx_t          mpoint )
{
  
}

//-----------------------------------------------------------------------
//  getHistoryNames
//-----------------------------------------------------------------------

void J2Plasticity::getHistoryNames   
    
  ( const StringVector&   hnames ) const

{

}


//-----------------------------------------------------------------------
//   Hist_ constructor ( empty )
//-----------------------------------------------------------------------

J2Plasticity::Hist_::Hist_()
{}