
/** @file J2Material.cpp
 *  @brief Implements a mixed isotropic-kinemetic linear hardening J2
 *  plasticity model.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 19 Aug 2024
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2024] 
 */

#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/util/Flex.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/LUSolver.h>

#include "util/Constants.h"
#include "util/Invariants.h"
#include "util/MathUtils.h"

#include "J2Material.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;
using jem::numeric::MatmulChain;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  J2Material::YOUNG_PROP          = "young";
const char*  J2Material::POISSON_PROP        = "poisson";
const char*  J2Material::YIELD_STRESS0_PROP  = "sig0";
const char*  J2Material::HARD_MODULUS_PROP   = "h";
const char*  J2Material::R_PROP              = "r";
const char*  J2Material::STATE_PROP          = "state";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


J2Material::J2Material 

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  young_   = 1.0;
  poisson_ = 0.0;
  sigY_    = 1.0;
  h_       = 1.0;
  r_       = 1.0;

  historyNames_.resize ( 9 );
  historyNames_[0] = "alf_xx";
  historyNames_[1] = "alf_yy";
  historyNames_[2] = "alf_zz";
  historyNames_[3] = "alf_xy";
  historyNames_[4] = "alf_yz";
  historyNames_[5] = "alf_xz";
  historyNames_[6] = "kap";
  historyNames_[7] = "mu";
  historyNames_[8] = "loading";
}


J2Material::~J2Material ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void J2Material::configure ( const Properties& props,
                             const Properties& globdat )
{
  props.get ( young_,   YOUNG_PROP         );
  props.get ( poisson_, POISSON_PROP       );
  props.get ( sigY_,    YIELD_STRESS0_PROP );
  props.get ( h_,       HARD_MODULUS_PROP  );
  props.get ( r_,       R_PROP, 0.0, 1.0   );

  G_ = young_ / 2.0 / ( 1.0 + poisson_ );
  K_ = young_ / 3.0 / ( 1.0 - 2.0 * poisson_ );

  // Compute strain count
  strCount_ = STRAIN_COUNTS[rank_];

  // Compute all projection matrices

  PVol_.resize ( strCount_, strCount_ ); PVol_ = 0.0;
  PDev_.resize ( strCount_, strCount_ ); PDev_ = 0.0;
  QDev_.resize ( strCount_, strCount_ ); QDev_ = 0.0;

  for ( int i = 0; i<3; i++)
  {
    for ( int j = 0; j<3; j++)
    {
      PVol_(i,j) = 1.0;
      i == j ? PDev_(i,j) = 2.0/3.0: PDev_(i,j) = -1.0/3.0;
      i == j ? QDev_(i,j) = 2.0/3.0: QDev_(i,j) = -1.0/3.0;
    }
  }

  if ( rank_ == 2 )
  {
    PDev_ (3,3) = 0.5;
    QDev_ (3,3) = 2.0;
  }
  else if ( rank_ == 3 )
  {
    PDev_(3,3)  = 0.5;
    PDev_(4,4)  = 0.5;
    PDev_(5,5)  = 0.5;

    QDev_   (3,3)  = 2.0;
    QDev_   (4,4)  = 2.0;
    QDev_   (5,5)  = 2.0;
  }

  // Compute the elastic stiffness matrix

  stiff0_.resize ( strCount_, strCount_ );
  stiff0_ = K_ * PVol_ + 2.0 * G_ * PDev_;
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void J2Material::getConfig ( const Properties& conf, 
                             const Properties& globdat ) const
{
  conf.set ( YOUNG_PROP        , young_  );
  conf.set ( POISSON_PROP      , poisson_);
  conf.set ( YIELD_STRESS0_PROP, sigY_   );
  conf.set ( HARD_MODULUS_PROP , h_      );
  conf.set ( R_PROP            , r_      );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void J2Material::updateConfig ()

{
  // If young_ and poisson_ are updated with setters, they must be updated
  // in elasticMat_, followed by a call to compute the updated material
  // stiffness matrix.

  G_ = young_ / 2.0 / ( 1.0 + poisson_ );
  K_ = young_ / 3.0 / ( 1.0 - 2.0 * poisson_ );

  // Compute the elastic stiffness matrix
  stiff0_ = K_ * PVol_ + 2.0 * G_ * PDev_;
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void J2Material::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )

{
  // Get the old step values (stress, strain, back-stress, kappa)
  Vector sig0 ( strCount_ );   sig0 = preHist_[ipoint].sig;
  Vector eps0 ( strCount_ );   eps0 = preHist_[ipoint].eps;
  Vector alf0 ( strCount_ );   alf0 = preHist_[ipoint].alf;
  double kap0;                 kap0 = preHist_[ipoint].kap;

  // Compute the strain increment
  Vector deps ( strCount_ );   deps = strain - eps0;

  // Compute the trial elastic stress
  Vector sigTr ( strCount_ ); sigTr = preHist_[ipoint].sig;
  sigTr += matmul( stiff0_, deps );

  // Compute the reduced trial elastic stress
  Vector sigX ( strCount_ ); 
  sigX = sigTr - alf0;

  // Compute the J2 stress
  const double sigE = sqrt( std::abs( 3.0/2.0*dot( sigX,matmul( QDev_,sigX ))));

  // Compute the yield function
  const double phi  = sigE - ( sigY_ + kap0 );

  // Check yield surface admissibility
  if ( phi < 0.0 )
  {
    // Elastic step: accept the trial solution
    stress = sigTr;
    stiff  = stiff0_;

    // Update the internal variables
    newHist_[ipoint].sig     = stress;
    newHist_[ipoint].eps     = strain;
    newHist_[ipoint].alf     = alf0;
    newHist_[ipoint].kap     = kap0;
    newHist_[ipoint].mu      = 0.0;
    newHist_[ipoint].loading = false;
  }
  else
  {
    // Plastic step: explicit solution

    // Some basic stuff
    MatmulChain<double,1>     mc1;

    // Compute the integrated plastic multiplier
    const double mu = phi / ( 3.0 * G_ + h_ );

    // Compute the flow direction (associative)
    Vector f ( strCount_ );    
    f   = (1.5 / sigE ) * mc1.matmul( QDev_,sigX );

    // Compute stress
    stress = sigTr - mu * mc1.matmul( stiff0_, f );

    // Compute material tangent stiffness
    stiff  = stiff0_ - 2.0 * G_ / sigE * ( 
             (2.0 * G_)/(3.0 * G_ + h_) * (kap0+sigY_)
             * matmul ( f, f )
             + 3.0 * mu *  G_ * PDev_
             );

    // Update internal variables
    newHist_[ipoint].sig     = stress;
    newHist_[ipoint].eps     = strain;
    newHist_[ipoint].alf     = alf0 + TWO_THIRD * ( 1.0 - r_) * f;
    newHist_[ipoint].kap     = kap0 + r_*h_*mu;
    newHist_[ipoint].mu      = mu;
    newHist_[ipoint].loading = true;
  }
}


//-----------------------------------------------------------------------
//  allocPoints
//-----------------------------------------------------------------------

void J2Material::allocPoints
    
  ( const idx_t          npoints )

{
  // Initialize history class

  Hist_ hist = Hist_();
  
  hist.sig .resize ( STRAIN_COUNTS[rank_] );
  hist.eps .resize ( STRAIN_COUNTS[rank_] );
  hist.alf .resize ( STRAIN_COUNTS[rank_] );

  hist.sig      = 0.0;
  hist.eps      = 0.0;
  hist.alf      = 0.0;
  hist.kap      = 0.0;
  hist.mu       = 0.0;
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

Ref<Material> J2Material::clone ( ) const

{
  return newInstance<J2Material> ( *this );
}


//-----------------------------------------------------------------------
//  commit
//-----------------------------------------------------------------------


void J2Material::commit ( )

{
  ( newHist_ ).swap ( preHist_ );

  latestHist_ = &( preHist_ );
  
}


//-----------------------------------------------------------------------
//  getHistory
//-----------------------------------------------------------------------


void J2Material::getHistory

  ( Vector&              hvals,
    const idx_t          ipoint )
{
  (*latestHist_)[ipoint].toVector ( hvals );
}

//-----------------------------------------------------------------------
//   Hist_ constructor ( empty )
//-----------------------------------------------------------------------

J2Material::Hist_::Hist_()
{}


// -------------------------------------------------------------------
//  Hist_::toVector
// -------------------------------------------------------------------

inline void J2Material::Hist_::toVector

 ( const Vector&  vec ) const

{
  vec[0] = alf[0];
  vec[1] = alf[1];
  vec[2] = alf[2];
  vec[3] = alf[3];
  vec[6] = kap;
  vec[7] = mu;
  vec[8] = int(loading);

  if ( alf.size() == 6 )
  {
    vec[4] = alf[4];
    vec[5] = alf[5];
  }
  else
  {
    vec[4] = 0.0;
    vec[5] = 0.0;
  }
}


// -------------------------------------------------------------------
//  Hist_::toPrint
// -------------------------------------------------------------------

inline void J2Material::Hist_::toPrint

 ( ) const

{
  System::out() <<   "alf " << alf
                << ", kap " << kap
                << ", mu "  << mu
                << "\n";
}