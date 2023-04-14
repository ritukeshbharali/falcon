#include <jem/io/PrintWriter.h>
#include <jem/base/Error.h>
#include <jem/base/System.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/utilities.h>

#include "HookeMaterial.h"
#include "util/BasicUtils.h"

#include <cstdlib>

using namespace jem;
using namespace jem::io;
using jem::numeric::matmul;


//=======================================================================
//   class HookeMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  HookeMaterial::YOUNG_PROP    = "young";
const char*  HookeMaterial::POISSON_PROP  = "poisson";
const char*  HookeMaterial::STATE_PROP    = "state";
const char*  HookeMaterial::AREA_PROP     = "area";

//-----------------------------------------------------------------------
//   constructor and destructor
//-----------------------------------------------------------------------

HookeMaterial::HookeMaterial

  ( idx_t       rank, const Properties& globdat )
    : Material ( rank, globdat )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  young_   = 1.0;
  poisson_ = 0.0;
  area_    = 1.0;

  stiffMat_.resize ( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  stiffMat_ = 0.0;
}

HookeMaterial::~HookeMaterial ()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void HookeMaterial::configure

  ( const Properties& props,
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

  if ( rank_ == 1 )
  {
    props.get ( area_, AREA_PROP );
  }

  // Get Problem type if material rank == 2

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

  computeStiffMat_ ();
}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void HookeMaterial::getConfig

  ( const Properties& conf,
    const Properties& globdat ) const
{
  conf.set ( YOUNG_PROP,   young_   );
  conf.set ( POISSON_PROP, poisson_ );
  conf.set ( STATE_PROP,   state_   );
  conf.set ( AREA_PROP,    area_    );
}

//-----------------------------------------------------------------------
//  allocPoints
//-----------------------------------------------------------------------

void HookeMaterial::allocPoints

  ( const idx_t       npoints )

{
}

//-----------------------------------------------------------------------
//  updateConfig
//-----------------------------------------------------------------------

void HookeMaterial::updateConfig ( )

{
  // updateConfig is run after the private/protected member is modified.
  // In this case, we recompute the stiffness matrix.

  computeStiffMat_ ();
}

//-----------------------------------------------------------------------
//  update
//-----------------------------------------------------------------------

void HookeMaterial::update

  ( Vector&       stress,
    Matrix&       stiff,
    const Vector& strain,
    const idx_t   ipoint )
{
  // get the elastic moduli

  stiff = stiffMat_;

  // compute stress

  matmul ( stress, stiff, strain );
}


//-----------------------------------------------------------------------
//   getStiffMat
//-----------------------------------------------------------------------


Matrix HookeMaterial::getStiffMat() const
{
  return stiffMat_;
}

//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> HookeMaterial::clone ( ) const
{
  return newInstance<HookeMaterial> ( *this );
}

//-----------------------------------------------------------------------
//  computeStiffMat_
//-----------------------------------------------------------------------

void   HookeMaterial::computeStiffMat_ () 
{
  const int     n  = STRAIN_COUNTS[rank_];

  const double  e  = young_;
  const double  nu = poisson_;


  if      ( rank_ == 1 )
  {
    stiffMat_(0,0) = e * area_;
  }
  else if ( rank_ == 3 )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(0,2) = a * nu;
    stiffMat_(1,1) = stiffMat_(0,0);
    stiffMat_(1,2) = stiffMat_(0,1);
    stiffMat_(2,2) = stiffMat_(0,0);
    stiffMat_(3,3) = a * b;
    stiffMat_(4,4) = stiffMat_(3,3);
    stiffMat_(5,5) = stiffMat_(3,3);

    // Copy lower triangle of the stress-strain matrix.

    for ( int i = 0; i < n; i++ )
    {
      for ( int j = 0; j < i; j++ )
      {
      stiffMat_(i,j) = stiffMat_(j,i);
      }
    }
  }
  else if ( state_ == PlaneStrain )
  {
    const double  a = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double  b = 0.5 * (1.0 - 2.0 * nu);
    const double  c = 1.0 - nu;

    stiffMat_(0,0) = a * c;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(0,2) = a * nu;
    stiffMat_(1,1) = a * c;
    stiffMat_(1,2) = a * nu;
    stiffMat_(3,3) = a * b;
    stiffMat_(2,0) = a * nu; stiffMat_(2,1) = a * nu; stiffMat_(2,2) = a * c;
  }
  else if ( state_ == PlaneStress )
  {
    const double  a = e / (1.0 - nu * nu);

    stiffMat_(0,0) = a;
    stiffMat_(0,1) = stiffMat_(1,0) = a * nu;
    stiffMat_(1,1) = a;
    stiffMat_(3,3) = a * 0.5 * (1.0 - nu);
  }
  else
  {
    throw Error ( JEM_FUNC, "unexpected rank: " + String ( rank_ ) );
  }

}