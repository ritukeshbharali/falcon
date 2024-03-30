#include <jem/base/Error.h>
#include <jem/util/Properties.h>

#include "Material.h"
#include "HookeMaterial.h"
#include "DamageMaterial.h"
#include "AmorPhaseMaterial.h"
#include "BourdinPhaseMaterial.h"
#include "MiehePhaseMaterial.h"

using namespace jem;


//=======================================================================
//   class Material
//=======================================================================

//-----------------------------------------------------------------------
//   constructors and destructor
//-----------------------------------------------------------------------

Material::Material

  ( const idx_t        rank,
    const Properties&  globdat )
{
  desperateMode_ = false;
  rank_          = rank;
}

Material::~Material()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void Material::configure

  ( const Properties& props,
    const Properties& globdat )
{}

//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------

void Material::getConfig

  ( const Properties& props,
    const Properties& globdat ) const
{}

//-----------------------------------------------------------------------
//   getHistory
//-----------------------------------------------------------------------

void Material::getHistory

  ( Vector&           hvals,
    const idx_t       mpoint )

{
  // Default implementation
  hvals.resize ( 0 );
}

//-----------------------------------------------------------------------
//   setHistory
//-----------------------------------------------------------------------

void Material::setHistory

  ( const Vector&    hvals,
    const idx_t      mpoint )

{
  // Default implementation
}

//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void Material::allocPoints

  ( const idx_t       npoints )

{}

//-----------------------------------------------------------------------
//   deallocPoints
//-----------------------------------------------------------------------

void Material::deallocPoints

  ( const idx_t       npoints )

{}

//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------

void Material::commit ()
{}

//-----------------------------------------------------------------------
//   checkCommit
//-----------------------------------------------------------------------

void Material::checkCommit

  ( const Properties& params )

{}

//-----------------------------------------------------------------------
//   commitOne
//-----------------------------------------------------------------------

void Material::commitOne

 ( const idx_t ipoint )

{}

//-----------------------------------------------------------------------
//   cancel
//-----------------------------------------------------------------------

void Material::cancel ()
{}

//-----------------------------------------------------------------------
//   pointCount
//-----------------------------------------------------------------------

idx_t  Material::pointCount () const

{
  return 0;
}

//-----------------------------------------------------------------------
//   isLoading
//-----------------------------------------------------------------------

bool  Material::isLoading

  ( idx_t ipoint ) const

{
  return false;
}

//-----------------------------------------------------------------------
//   wasLoading
//-----------------------------------------------------------------------

bool  Material::wasLoading

  ( idx_t ipoint ) const

{
  return false;
}

//-----------------------------------------------------------------------
//   despair
//-----------------------------------------------------------------------

bool Material::despair ()

{
  desperateMode_ = true;

  idx_t np = pointCount();

  hasSwitched_.resize ( np );
  useSecant_  .resize ( np );

  useSecant_ = false;

  for ( idx_t ip = 0; ip < np; ++ip )
  {
    hasSwitched_[ip] = ( isLoading(ip) != wasLoading(ip) );
  }

  return np > 0;
}

//-----------------------------------------------------------------------
//   endDespair
//-----------------------------------------------------------------------

void Material::endDespair ()

{
  desperateMode_ = false;
  useSecant_     = false;
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newMaterial
//-----------------------------------------------------------------------

Ref<Material>  newMaterial

  ( const String&     name,
    const Properties& conf,
    const Properties& props,
    const Properties& globdat )

{
  Properties    matProps = props.getProps ( name );
  Properties    matConf  = conf.makeProps ( name );

  Ref<Material> mat;
  String        type;
  idx_t         rank;

  matProps.get ( type, "type" );
  matConf.set  ( "type", type );

  matProps.get ( rank, "rank" );
  matConf.set  ( "rank", rank );

  if      ( type == "Hooke" )
    mat = newInstance<HookeMaterial>        ( rank, globdat );
  else if ( type == "Damage" )
    mat = newInstance<DamageMaterial>       ( rank, globdat );
  else if ( type == "AmorPhase" )
    mat = newInstance<AmorPhaseMaterial>    ( rank, globdat );
  else if ( type == "BourdinPhase" )
    mat = newInstance<BourdinPhaseMaterial> ( rank, globdat );
  else if ( type == "MiehePhase" )
    mat = newInstance<MiehePhaseMaterial>   ( rank, globdat );
  else
    matProps.propertyError ( name, "Invalid material: " + type );

  return mat;
}