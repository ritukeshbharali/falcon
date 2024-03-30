
/** @file Invariants.cpp
 *  @brief Computes the invariants.
 * 
 *  @note: Output is always a Tuple vector.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */

/* Include jem and jive headers */

#include <jem/base/Array.h>
#include <jive/util/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/EigenUtils.h>

/* Include other headers */

#include "Constants.h"
#include "Invariants.h"
#include "MathUtils.h"

using namespace std;
using namespace jem;

//=======================================================================
//   class Invariants
//=======================================================================

//-----------------------------------------------------------------------------
//  constructor and destructor
//-----------------------------------------------------------------------------

Invariants::Invariants 

    ( const Vector& vec,
      const bool    isStrain )
{
  JEM_ASSERT ( vec.size() == 4 || vec.size() == 6 );

  // Factor 2 is removed when dealing with strains.

  double fact = ( isStrain ) ? 0.5 : 1.0;

  // 2D (four components)

  v6_[0] = vec[0];
  v6_[1] = vec[1];
  v6_[2] = vec[2];
  v6_[3] = fact * vec[3];
  v6_[4] = 0.0;
  v6_[5] = 0.0;

  // 3D (six components)

  if ( vec.size() == 6 )
  {
    v6_[4] = fact * vec[4];
    v6_[5] = fact * vec[5];
  }
}


//-----------------------------------------------------------------------------
//  getI1
//-----------------------------------------------------------------------------

double Invariants::getI1 ()

{  
  double I1 = v6_[0] + v6_[1] + v6_[2];
  
  return I1;  
}

//-----------------------------------------------------------------------------
//  getI2
//-----------------------------------------------------------------------------

double Invariants::getI2 ()

{  
  double I2 = v6_[0] * v6_[1] + v6_[1] * v6_[2] + v6_[2] * v6_[0] - 
            ( v6_[3] * v6_[3] + v6_[4] * v6_[4] + v6_[5] * v6_[5] );

  return I2;  
}


//-----------------------------------------------------------------------------
//  getI3
//-----------------------------------------------------------------------------

double Invariants::getI3 ()

{
  double I3 = v6_[0] * v6_[1] * v6_[2] + 2.0 * v6_[3] * v6_[4] * v6_[5] -
              v6_[0] * v6_[4] * v6_[4] - 
              v6_[1] * v6_[5] * v6_[5] - 
              v6_[2] * v6_[3] * v6_[3];

  return I3;  
}


//-----------------------------------------------------------------------------
//  getJ2
//-----------------------------------------------------------------------------

double Invariants::getJ2 ()

{
  double J2 = ONE_THIRD * ( v6_[0] * v6_[0] + v6_[1] * v6_[1] + v6_[2] * v6_[2]   -
                     v6_[0] * v6_[1] - v6_[1] * v6_[2] - v6_[2] * v6_[0] ) +
                     v6_[3] * v6_[3] + v6_[4] * v6_[4] + v6_[5] * v6_[5];
  return J2; 
}


// ---------------------------------------------------------------------
//  getGradI1
// ---------------------------------------------------------------------

Vec6 Invariants::getGradI1 () 
  
{
  Vec6 gradI1;

  gradI1[0] = 1.0;
  gradI1[1] = 1.0;
  gradI1[2] = 1.0;
  
  gradI1[3] = 0.0;
  gradI1[4] = 0.0;
  gradI1[5] = 0.0;    

  return gradI1;
}


// ---------------------------------------------------------------------
//  getGradI2
// ---------------------------------------------------------------------

Vec6  Invariants::getGradI2 ()

{
  Vec6 gradI2;

  gradI2[0] = v6_[1] + v6_[2];
  gradI2[1] = v6_[2] + v6_[0];
  gradI2[2] = v6_[0] + v6_[1];

  gradI2[3] = - 2. * v6_[3];
  gradI2[4] = - 2. * v6_[4];
  gradI2[5] = - 2. * v6_[5];

  return gradI2;
}


// ---------------------------------------------------------------------
//  getGradI3
// ---------------------------------------------------------------------

Vec6  Invariants::getGradI3 ()
  
{
  Vec6 gradI3;

  gradI3[0] = v6_[1] * v6_[2] - v6_[4] * v6_[4];
  gradI3[1] = v6_[2] * v6_[0] - v6_[5] * v6_[5];
  gradI3[2] = v6_[0] * v6_[1] - v6_[3] * v6_[3];

  gradI3[3] = 2.0 * ( v6_[4] * v6_[5] - v6_[2] * v6_[3] );
  gradI3[4] = 2.0 * ( v6_[5] * v6_[3] - v6_[0] * v6_[4] );
  gradI3[5] = 2.0 * ( v6_[3] * v6_[4] - v6_[1] * v6_[5] );

  return gradI3;
}


// ---------------------------------------------------------------------
//  getGradJ2
// ---------------------------------------------------------------------

Vec6  Invariants::getGradJ2 ()

{
  Vec6 gradJ2;

  gradJ2[0] = ONE_THIRD * ( 2.0 * v6_[0] - v6_[1] - v6_[2] );
  gradJ2[1] = ONE_THIRD * ( 2.0 * v6_[1] - v6_[2] - v6_[0] );
  gradJ2[2] = ONE_THIRD * ( 2.0 * v6_[2] - v6_[0] - v6_[1] );

  gradJ2[3] = 2.0 * v6_[3];
  gradJ2[4] = 2.0 * v6_[4];
  gradJ2[5] = 2.0 * v6_[5];

  return gradJ2;
}


// ---------------------------------------------------------------------
//  getPrincipalValues
// ---------------------------------------------------------------------

Vec3  Invariants::getPrincipalValues ()

{
  const double I1 = getI1 ();
  const double I2 = getI2 ();
  const double I3 = getI3 ();

  Vec3 roots;

  const idx_t  n  = solveCubicEqn ( roots, 1.0, -I1, I2, -I3 );

  for ( idx_t i = n; i < 3; ++i ) 
    roots[i] = 0.0;

  return roots;
}


void Invariants::getPrincipalValues 
  
  ( const Vector& princ )

{
  JEM_ASSERT ( princ.size() == 3 );

  Vec3 roots = getPrincipalValues ();

  princ[0] = roots[0];
  princ[1] = roots[1];
  princ[2] = roots[2];
}