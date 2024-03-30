
/** @file TupleUtils.h
 *  @brief Implements operations on Tuples (fixed size arrays).
 *         Originally written by Frans van der Meer.
 *  
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */

#ifndef TUPLE_UTILS_H
#define TUPLE_UTILS_H

#include "Arrays.h"

//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------

void            m6_to_tt6

  (       Mat6&   tt6,
    const Matrix& m6  );
      
void            tt6_to_m6

  ( const Matrix& m6,
    const Mat6&   tt6 );
      
void v6_to_t6

  (       Vec6&   t6,
    const Vector& v6 );
    
Vec6 v6_to_t6

  ( const Vector& v6 );
  
void t6_to_v6

  ( const Vector& v6,
    const Vec6&   t6 );

void t6_to_v4

  ( const Vector& v4,
    const Vec6&   t6 );     
   
Vec6 tmatmul

  ( const Matrix& mat,
    const Vec6&   t6 );
    
Vec6 tmatmul

  ( const Vec6&   t6,
    const Matrix& mat );
    
Vec6 tmatmul

  ( const Vec6&   t6,
    const Mat6&   tt6 );
    
Vec6 tmatmul

  ( const Mat6&   tt6,
    const Vec6&   t6 );
    
Mat6 tmatmul

  ( const Vec6&   t0,
    const Vec6&   t1  );
    
Mat6 tmatmul

  ( const Mat6&     tt6,
    const Matrix&   mat  );    
    
Mat6 tmatmul

  ( const Matrix&   mat,
    const Mat6&     tt6  ); 
    
Mat6 tmatmul

  ( const Mat6&     tt6_1,
    const Mat6&     tt6_2  ); 
    
double tdot

  ( const Vector&   v6, 
    const Vec6&     t6 );    

double tdot

  ( const Vec6&     t6, 
    const Vector&   v6 );

#endif
