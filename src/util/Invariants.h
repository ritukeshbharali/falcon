
/** @file Invariants.h
 *  @brief Computes the invariants.
 * 
 *  @note: Output is always a Tuple vector, size 6.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */

#ifndef INVARIANTS_H
#define INVARIANTS_H

/* Include jem and jive headers */

#include "Arrays.h"

//=======================================================================
//   class Invariants
//=======================================================================

/** @brief 
 *  The Invariants class computes the invariants of the stress and 
 *  strain.
 */

class Invariants
{
 public:

  Invariants ( const Vector& vec,
               const bool    isStrain );

  // Compute invariants

  double getI1     ();
  double getI2     ();
  double getI3     ();
  double getJ2     ();

  // Compute derivatives of invariants

  Vec6 getGradI1 ();
  Vec6 getGradI2 ();
  Vec6 getGradI3 ();
  Vec6 getGradJ2 ();

  // Compute the principal values

  Vec3   getPrincipalValues ();
  void   getPrincipalValues ( const Vector& princ );

 private:

  Vec6  v6_;
  
};

#endif
