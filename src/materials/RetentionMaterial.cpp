#include <jem/base/Error.h>
#include <jem/util/Properties.h>

#include "RetentionMaterial.h"
#include "VanGenuchtenRetentionMaterial.h"
#include "LiakopoulosRetentionMaterial.h"

using namespace jem;


//=======================================================================
//   class RetentionMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   constructors and destructor
//-----------------------------------------------------------------------

RetentionMaterial::RetentionMaterial

  ( const Properties&  globdat )
{}

RetentionMaterial::~RetentionMaterial()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void RetentionMaterial::configure

  ( const Properties& props,
    const Properties& globdat )
{}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void RetentionMaterial::getConfig

  ( const Properties& props,
    const Properties& globdat ) const
{}

//-----------------------------------------------------------------------
//   update (default implementation)
//-----------------------------------------------------------------------

void RetentionMaterial::updateSaturation

    ( double&              Sf,
      double&              dSfdpf,
      const double&        pf      )
{}    

void RetentionMaterial::updateSaturation

    ( double&              Sf,
      double&              dSfdpf,
      const double&        pf,
      const double&        pg      )      
{}    

void RetentionMaterial::updatePermeability

    ( double&              krf,
      double&              dkrfdpf,
      const double&        pf      )
{}    

void RetentionMaterial::updatePermeability

    ( double&              krf,
      double&              dkrfdpf,
      const double&        pf,
      const double&        pg      )        
{}    

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newRetentionMaterial
//-----------------------------------------------------------------------

Ref<RetentionMaterial>  newRetentionMaterial

  ( const String&     name,
    const Properties& conf,
    const Properties& props,
    const Properties& globdat )

{
  Properties    matProps = props.getProps ( name );
  Properties    matConf  = conf.makeProps ( name );

  Ref<RetentionMaterial> mat;
  String        type;

  matProps.get ( type, "type" );
  matConf.set  ( "type", type );

  if ( type == "VanGenuchten" )
  {
    mat = newInstance<VanGenuchtenRetentionMaterial> ( globdat );
  }
  else if ( type == "Liakopoulos" )
  {
    mat = newInstance<LiakopoulosRetentionMaterial> ( globdat );
  }
  else
  {
    matProps.propertyError (
      name,
      "invalid retention material type: " + type
    );
  }

  return mat;
}