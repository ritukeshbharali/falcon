
/** @file LiakopoulosRetentionMaterial.cpp
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

#include "util/BasicUtils.h"
#include "LiakopoulosRetentionMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  LiakopoulosRetentionMaterial::RHO_PROP  = "rho";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


LiakopoulosRetentionMaterial::LiakopoulosRetentionMaterial 

  ( const Properties& globdat )
    : RetentionMaterial ( globdat )
{
  rho_  = 1000.0;
  g_    = 9.81;
}


LiakopoulosRetentionMaterial::~LiakopoulosRetentionMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void LiakopoulosRetentionMaterial::configure ( const Properties& props,
                                const Properties& globdat )
{
  props.find ( rho_,  RHO_PROP  );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void LiakopoulosRetentionMaterial::getConfig ( const Properties& conf, 
                                const Properties& globdat ) const
{
  conf.set ( RHO_PROP , rho_  );
}

//-----------------------------------------------------------------------
//   update (single fluid phase)
//-----------------------------------------------------------------------

void LiakopoulosRetentionMaterial::update

    ( double&               Sf,
      double&               dSfdpf,
      double&               krf,
      double&               dkrfdpf,
      const double&         pf       )

{
  const double H     = -pf / ( rho_ * g_  );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /**
     * compute saturation and its derivatives 
     * DOI: 10.3970/cmes.2009.053.255 (Eqn.30)
     */
               
    Sf      = 1.0 - 1.9722e-11 * ( ::pow( -pf, 2.4279 ) );
    dSfdpf  =    4.7883044e-11 * ( ::pow( -pf, 1.4279 ) );

    krf     = 1.0 - 2.207 * ( ::pow( 1.0-Sf, 0.9529  ) );
    dkrfdpf =   2.1030503 * ( ::pow( 1.0-Sf, -0.0471 ) );

 }
}

//-----------------------------------------------------------------------
//   update (fluid and gas phases present)
//-----------------------------------------------------------------------

void LiakopoulosRetentionMaterial::update

    ( double&               Sf,
      double&               dSfdpf,
      double&               krf,
      double&               dkrfdpf,
      const double&         pf,
      const double&         pg      )

{
  const double H     = ( pg - pf ) / ( rho_ * g_  );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /**
     * compute saturation and its derivatives 
     * DOI: 10.3970/cmes.2009.053.255 (Eqn.30)
     */
               
    Sf      = 1.0 - 1.9722e-11 * ( ::pow( (pg-pf), 2.4279 ) );
    dSfdpf  =    4.7883044e-11 * ( ::pow( (pg-pf), 1.4279 ) );

    krf     = 1.0 - 2.207 * ( ::pow( 1.0-Sf, 0.9529  ) );
    dkrfdpf =   2.1030503 * ( ::pow( 1.0-Sf, -0.0471 ) );
 }
}

//-----------------------------------------------------------------------
//   updateSaturation (single fluid phase)
//-----------------------------------------------------------------------

void LiakopoulosRetentionMaterial::updateSaturation

    ( double&               Sf,
      double&               dSfdpf,
      const double&         pf       )

{
  const double H     = -pf / ( rho_ * g_  );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
  }
  else
  {
    /**
     * compute saturation and its derivatives 
     * DOI: 10.3970/cmes.2009.053.255 (Eqn.30)
     */
               
    Sf      = 1.0 - 1.9722e-11 * ( ::pow( -pf, 2.4279 ) );
    dSfdpf  =    4.7883044e-11 * ( ::pow( -pf, 1.4279 ) );
 }
}

//-----------------------------------------------------------------------
//   updateSaturation (fluid and gas phases present)
//-----------------------------------------------------------------------

void LiakopoulosRetentionMaterial::updateSaturation

    ( double&               Sf,
      double&               dSfdpf,
      const double&         pf,
      const double&         pg      )

{
  const double H     = ( pg - pf ) / ( rho_ * g_  );

  if ( H < 1e-20 )
  {
    Sf       = 1.0;
    dSfdpf   = 0.0;
  }
  else
  {
    /**
     * compute saturation and its derivatives 
     * DOI: 10.3970/cmes.2009.053.255 (Eqn.30)
     */
               
    Sf      = 1.0 - 1.9722e-11 * ( ::pow( (pg-pf), 2.4279 ) );
    dSfdpf  =    4.7883044e-11 * ( ::pow( (pg-pf), 1.4279 ) );
 }
}

//-----------------------------------------------------------------------
//   updatePermeability (single fluid phase)
//-----------------------------------------------------------------------

void LiakopoulosRetentionMaterial::updatePermeability

    ( double&               krf,
      double&               dkrfdpf,
      const double&         pf       )

{
  const double H     = -pf / ( rho_ * g_  );

  if ( H < 1e-20 )
  {
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /**
     * compute saturation and its derivatives 
     * DOI: 10.3970/cmes.2009.053.255 (Eqn.30)
     */
               
    double Sf = 1.0 - 1.9722e-11 * ( ::pow( -pf, 2.4279 ) );

    krf       = 1.0 - 2.207 * ( ::pow( 1.0-Sf, 0.9529  ) );
    dkrfdpf   =   2.1030503 * ( ::pow( 1.0-Sf, -0.0471 ) );

 }
}

//-----------------------------------------------------------------------
//   updatePermeability (fluid and gas phases present)
//-----------------------------------------------------------------------

void LiakopoulosRetentionMaterial::updatePermeability

    ( double&               krf,
      double&               dkrfdpf,
      const double&         pf,
      const double&         pg      )

{
  const double H     = ( pg - pf ) / ( rho_ * g_  );

  if ( H < 1e-20 )
  {
    krf      = 1.0;
    dkrfdpf  = 0.0;
  }
  else
  {
    /**
     * compute saturation and its derivatives 
     * DOI: 10.3970/cmes.2009.053.255 (Eqn.30)
     */
               
    double Sf = 1.0 - 1.9722e-11 * ( ::pow( (pg-pf), 2.4279 ) );

    krf       = 1.0 - 2.207 * ( ::pow( 1.0-Sf, 0.9529  ) );
    dkrfdpf   =   2.1030503 * ( ::pow( 1.0-Sf, -0.0471 ) );
 }
}

//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<RetentionMaterial> LiakopoulosRetentionMaterial::clone ( ) const

{
  return newInstance<LiakopoulosRetentionMaterial> ( *this );
}