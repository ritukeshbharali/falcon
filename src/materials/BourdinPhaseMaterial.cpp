
/** @file BourdinPhaseMaterial.cpp
 *  @brief Implements a phase-field fracture material model with no split.
 *  
 *  This class implements a phase-field fracture material
 *  model without any split in the strain energy density
 *  (see DOI: 10.1016/S0022-5096(99)00028-9).
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
#include "BourdinPhaseMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  BourdinPhaseMaterial::YOUNG_PROP          = "young";
const char*  BourdinPhaseMaterial::POISSON_PROP        = "poisson";
const char*  BourdinPhaseMaterial::STATE_PROP          = "state";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


BourdinPhaseMaterial::BourdinPhaseMaterial 

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  elasticMat_ = newInstance<HookeMaterial> ( rank, globdat );
  
  elasticStiffMat_.resize( STRAIN_COUNTS[rank_], STRAIN_COUNTS[rank_] );
  elasticStiffMat_ = 0.0;

  stressP_.resize( STRAIN_COUNTS[rank_] );
  stressP_ = 0.0;
}


BourdinPhaseMaterial::~BourdinPhaseMaterial ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void BourdinPhaseMaterial::configure ( const Properties& props,
                                const Properties& globdat )
{
  elasticMat_->configure ( props, globdat );
  
  elasticStiffMat_ = elasticMat_->getStiffMat();
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void BourdinPhaseMaterial::getConfig ( const Properties& conf, 
                                const Properties& globdat ) const
{
  elasticMat_->getConfig ( conf, globdat );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void BourdinPhaseMaterial::updateConfig ()

{
  // If young_ and poisson_ are updated with setters, they must be updated
  // in elasticMat_, followed by a call to compute the updated material
  // stiffness matrix.
}


//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void BourdinPhaseMaterial::update

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

void BourdinPhaseMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const double&         gphi,
      int                   ipoint )

{
  stressP_ = matmul( elasticStiffMat_,strain );
  Psi_     = 0.5 * dot( stressP_,strain );
  
  stress   = gphi * stressP_;
  stiff    = gphi * elasticStiffMat_;
}


//-----------------------------------------------------------------------
//   update (overloaded for micromorphic variants)
//-----------------------------------------------------------------------

void BourdinPhaseMaterial::update

    ( Vector&               stressP,
      Vector&               stressN,
      Matrix&               stiffP,
      Matrix&               stiffN,
      const Vector&         strain,
      int                   ip     )

{

  stressP  = matmul( elasticStiffMat_,strain );
  stressN  = 0.0;

  stiffP   = elasticStiffMat_;
  stiffN   = 0.0;

  stressP_ = stressP;
  Psi_     = 0.5 * dot( stressP_,strain );

}


//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> BourdinPhaseMaterial::clone ( ) const

{
  return newInstance<BourdinPhaseMaterial> ( *this );
}