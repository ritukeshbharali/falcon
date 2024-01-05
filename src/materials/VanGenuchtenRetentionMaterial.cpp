
/** @file VanGenuchtenRetentionMaterial.cpp
 *  @brief Implements a van Genuchten fluid retention model.
 *  
 *  This class implements van Genuchten fluid retention model. 
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 20 October 2022
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2022]
 */


#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/utilities.h>

#include "util/Constants.h"
#include "VanGenuchtenRetentionMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  VanGenuchtenRetentionMaterial::A_PROP    = "aVG";
const char*  VanGenuchtenRetentionMaterial::E_PROP    = "eVG";
const char*  VanGenuchtenRetentionMaterial::J_PROP    = "jVG";
const char*  VanGenuchtenRetentionMaterial::SRES_PROP = "resSaturation";
const char*  VanGenuchtenRetentionMaterial::RHO_PROP  = "rho";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


VanGenuchtenRetentionMaterial::VanGenuchtenRetentionMaterial 

  ( const Properties& globdat )
    : RetentionMaterial ( globdat )
{
  aVG_  = 0.028;
  eVG_  = 0.5;
  jVG_  = 1.3;
  sRes_ = 0.15;
  rho_  = 1000.0;
  g_    = 9.81;
}


VanGenuchtenRetentionMaterial::~VanGenuchtenRetentionMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void VanGenuchtenRetentionMaterial::configure ( const Properties& props,
                                const Properties& globdat )
{
  props.find ( aVG_,  A_PROP    );
  props.find ( eVG_,  E_PROP    );
  props.find ( jVG_,  J_PROP    );
  props.find ( sRes_, SRES_PROP );
  props.find ( rho_,  RHO_PROP  );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void VanGenuchtenRetentionMaterial::getConfig ( const Properties& conf, 
                                const Properties& globdat ) const
{
  conf.set ( A_PROP   , aVG_  );
  conf.set ( E_PROP   , eVG_  );
  conf.set ( J_PROP   , jVG_  );
  conf.set ( SRES_PROP, sRes_ );
  conf.set ( RHO_PROP , rho_  );
}

//-----------------------------------------------------------------------
//   update (single fluid phase)
//-----------------------------------------------------------------------

void VanGenuchtenRetentionMaterial::update

    ( double&               Sf,
      double&               dSfdpf,
      double&               krf,
      double&               dkrfdpf,
      const double&         pf       )

{
  const double H     = -pf / ( rho_ * g_  );
  const double hVG_  = 1.0 - ( 1.0 / jVG_ );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /*
     * compute saturation and its derivatives using van Genuchten
     * relations, DOI: 10.2136/sssaj1980.03615995004400050002x
     */

    double tmp    = ::pow( aVG_ * H, jVG_ );
    double Se     = ::pow( 1.0 + tmp, -hVG_ );
               
    Sf            = ( 1.0 - sRes_ ) * Se + sRes_;

    double dSfdSe = ( 1.0 - sRes_ );
    double dSedH  = - jVG_ * hVG_ * ::pow( aVG_ * H, jVG_ - 1.0 ) * ::pow( 1.0 + tmp, -hVG_-1.0 ) * aVG_; 
    double dHdpf  = - 1.0 / ( rho_ * g_ );

    dSfdpf        = dSfdSe * dSedH * dHdpf;

    /*
     * compute relative permeability and its derivative using Mualem
     * relations, DOI: 10.1029/WR012i003p00513
     */

    krf = ::pow( Se, eVG_ ) * ::pow( ( 1.0 - ::pow( ( 1.0 - ::pow( Se, 1.0/hVG_ ) ), hVG_ ) ), 2.0 );
    tmp = ::pow( Se, 1./hVG_ );

    double dkrfdSe = (::pow( ::pow( 1. - tmp, hVG_ ) - 1. , 2. ) 
                       - 4. * tmp * ::pow( 1. - tmp, hVG_ - 1. ) 
                       * ( ::pow( 1. - tmp, hVG_ ) - 1. ) ) / 2. * ::sqrt(Se);

    dkrfdpf  = dkrfdSe * dSedH * dHdpf; 

 }
}

//-----------------------------------------------------------------------
//   update (fluid and gas phases present)
//-----------------------------------------------------------------------

void VanGenuchtenRetentionMaterial::update

    ( double&               Sf,
      double&               dSfdpf,
      double&               krf,
      double&               dkrfdpf,
      const double&         pf,
      const double&         pg      )

{
  const double H     = ( pg - pf ) / ( rho_ * g_  );
  const double hVG_  = 1.0 - ( 1.0 / jVG_ );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /*
     * compute saturation and its derivatives using van Genuchten
     * relations, DOI: 10.2136/sssaj1980.03615995004400050002x
     */

    double tmp    = ::pow( aVG_ * H, jVG_ );
    double Se     = ::pow( 1.0 + tmp, -hVG_ );
               
    Sf            = ( 1.0 - sRes_ ) * Se + sRes_;

    double dSfdSe = ( 1.0 - sRes_ );
    double dSedH  = - jVG_ * hVG_ * ::pow( aVG_ * H, jVG_ - 1.0 ) * ::pow( 1.0 + tmp, -hVG_-1.0 ) * aVG_; 
    double dHdpf  = - 1.0 / ( rho_ * g_ );

    dSfdpf        = dSfdSe * dSedH * dHdpf;

    /*
     * compute relative permeability and its derivative using Mualem
     * relations, DOI: 10.1029/WR012i003p00513
     */

    krf = ::pow( Se, eVG_ ) * ::pow( ( 1.0 - ::pow( ( 1.0 - ::pow( Se, 1.0/hVG_ ) ), hVG_ ) ), 2.0 );
    tmp = ::pow( Se, 1./hVG_ );

    double dkrfdSe = (::pow( ::pow( 1. - tmp, hVG_ ) - 1. , 2. ) 
                       - 4. * tmp * ::pow( 1. - tmp, hVG_ - 1. ) 
                       * ( ::pow( 1. - tmp, hVG_ ) - 1. ) ) / 2. * ::sqrt(Se);

    dkrfdpf  = dkrfdSe * dSedH * dHdpf; 
 }
}

//-----------------------------------------------------------------------
//   updateSaturation (single fluid phase)
//-----------------------------------------------------------------------

void VanGenuchtenRetentionMaterial::updateSaturation

    ( double&               Sf,
      double&               dSfdpf,
      const double&         pf       )

{
  const double H     = -pf / ( rho_ * g_  );
  const double hVG_  = 1.0 - ( 1.0 / jVG_ );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
  }
  else
  {
    /*
     * compute saturation and its derivatives using van Genuchten
     * relations, DOI: 10.2136/sssaj1980.03615995004400050002x
     */

    double tmp    = ::pow( aVG_ * H, jVG_ );
    double Se     = ::pow( 1.0 + tmp, -hVG_ );
               
    Sf            = ( 1.0 - sRes_ ) * Se + sRes_;

    double dSfdSe = ( 1.0 - sRes_ );
    double dSedH  = - jVG_ * hVG_ * ::pow( aVG_ * H, jVG_ - 1.0 ) * ::pow( 1.0 + tmp, -hVG_-1.0 ) * aVG_; 
    double dHdpf  = - 1.0 / ( rho_ * g_ );

    dSfdpf        = dSfdSe * dSedH * dHdpf;
 }
}

//-----------------------------------------------------------------------
//   updateSaturation (fluid and gas phases present)
//-----------------------------------------------------------------------

void VanGenuchtenRetentionMaterial::updateSaturation

    ( double&               Sf,
      double&               dSfdpf,
      const double&         pf,
      const double&         pg      )

{
  const double H     = ( pg - pf ) / ( rho_ * g_  );
  const double hVG_  = 1.0 - ( 1.0 / jVG_ );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
  }
  else
  {
    /*
     * compute saturation and its derivatives using van Genuchten
     * relations, DOI: 10.2136/sssaj1980.03615995004400050002x
     */

    double tmp    = ::pow( aVG_ * H, jVG_ );
    double Se     = ::pow( 1.0 + tmp, -hVG_ );
               
    Sf            = ( 1.0 - sRes_ ) * Se + sRes_;

    double dSfdSe = ( 1.0 - sRes_ );
    double dSedH  = - jVG_ * hVG_ * ::pow( aVG_ * H, jVG_ - 1.0 ) * ::pow( 1.0 + tmp, -hVG_-1.0 ) * aVG_; 
    double dHdpf  = - 1.0 / ( rho_ * g_ );

    dSfdpf        = dSfdSe * dSedH * dHdpf;
 }
}

//-----------------------------------------------------------------------
//   updatePermeability (single fluid phase)
//-----------------------------------------------------------------------

void VanGenuchtenRetentionMaterial::updatePermeability

    ( double&               krf,
      double&               dkrfdpf,
      const double&         pf       )

{
  const double H     = -pf / ( rho_ * g_  );
  const double hVG_  = 1.0 - ( 1.0 / jVG_ );

  if ( H < 1e-20 )
  {
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /*
     * compute saturation and its derivatives using van Genuchten
     * relations, DOI: 10.2136/sssaj1980.03615995004400050002x
     */

    double tmp    = ::pow( aVG_ * H, jVG_ );
    double Se     = ::pow( 1.0 + tmp, -hVG_ );
    double dSedH  = - jVG_ * hVG_ * ::pow( aVG_ * H, jVG_ - 1.0 ) * ::pow( 1.0 + tmp, -hVG_-1.0 ) * aVG_; 
    double dHdpf  = - 1.0 / ( rho_ * g_ );

    /*
     * compute relative permeability and its derivative using Mualem
     * relations, DOI: 10.1029/WR012i003p00513
     */

    krf = ::pow( Se, eVG_ ) * ::pow( ( 1.0 - ::pow( ( 1.0 - ::pow( Se, 1.0/hVG_ ) ), hVG_ ) ), 2.0 );
    tmp = ::pow( Se, 1./hVG_ );

    double dkrfdSe = (::pow( ::pow( 1. - tmp, hVG_ ) - 1. , 2. ) 
                       - 4. * tmp * ::pow( 1. - tmp, hVG_ - 1. ) 
                       * ( ::pow( 1. - tmp, hVG_ ) - 1. ) ) / 2. * ::sqrt(Se);

    dkrfdpf  = dkrfdSe * dSedH * dHdpf; 

 }
}

//-----------------------------------------------------------------------
//   updatePermeability (fluid and gas phases present)
//-----------------------------------------------------------------------

void VanGenuchtenRetentionMaterial::updatePermeability

    ( double&               krf,
      double&               dkrfdpf,
      const double&         pf,
      const double&         pg      )

{
  const double H     = ( pg - pf ) / ( rho_ * g_  );
  const double hVG_  = 1.0 - ( 1.0 / jVG_ );

  if ( H < 1e-20 )
  {
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /*
     * compute saturation and its derivatives using van Genuchten
     * relations, DOI: 10.2136/sssaj1980.03615995004400050002x
     */

    double tmp    = ::pow( aVG_ * H, jVG_ );
    double Se     = ::pow( 1.0 + tmp, -hVG_ );
    double dSedH  = - jVG_ * hVG_ * ::pow( aVG_ * H, jVG_ - 1.0 ) * ::pow( 1.0 + tmp, -hVG_-1.0 ) * aVG_; 
    double dHdpf  = - 1.0 / ( rho_ * g_ );

    /*
     * compute relative permeability and its derivative using Mualem
     * relations, DOI: 10.1029/WR012i003p00513
     */

    krf = ::pow( Se, eVG_ ) * ::pow( ( 1.0 - ::pow( ( 1.0 - ::pow( Se, 1.0/hVG_ ) ), hVG_ ) ), 2.0 );
    tmp = ::pow( Se, 1./hVG_ );

    double dkrfdSe = (::pow( ::pow( 1. - tmp, hVG_ ) - 1. , 2. ) 
                       - 4. * tmp * ::pow( 1. - tmp, hVG_ - 1. ) 
                       * ( ::pow( 1. - tmp, hVG_ ) - 1. ) ) / 2. * ::sqrt(Se);

    dkrfdpf  = dkrfdSe * dSedH * dHdpf; 
 }
}

//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<RetentionMaterial> VanGenuchtenRetentionMaterial::clone ( ) const

{
  return newInstance<VanGenuchtenRetentionMaterial> ( *this );
}