
/** @file TupleUtils.cpp
 *  @brief Implements operations on Tuples (fixed size arrays).
 *         Originally written by Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */


#include <jem/base/array/select.h>
#include <jem/base/array/utilities.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/Error.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/matmul.h>

#include "TupleUtils.h"

using jem::ALL;
using jem::numeric::matmul;


//-----------------------------------------------------------------------
//   functions - Array <-> Tuple conversions
//-----------------------------------------------------------------------

void m6_to_tt6

  (       Mat6&   tt6,
    const Matrix& m6   )

{
  for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        tt6(i,j) = m6(i,j);
      }
    };
}

void tt6_to_m6

  ( const Matrix& m6,
    const Mat6&   tt6 )

{
  for ( idx_t i = 0; i < 6; ++i )
    {
      for ( idx_t j = 0; j < 6; ++j )
      {
        m6(i,j) = tt6(i,j);
      }
    };
}

void v6_to_t6

  (       Vec6&   t6,
    const Vector& v6 )

{
  t6[0] = v6[0]; t6[1] = v6[1]; t6[2] = v6[2]; 
  t6[3] = v6[3]; t6[4] = v6[4]; t6[5] = v6[5]; 
}

Vec6 v6_to_t6

  ( const Vector& v6 )

{
  Vec6 t6;
  v6_to_t6 ( t6, v6 );
  return t6;
}

void t6_to_v6

  ( const Vector& v6,
    const Vec6&   t6 )

{
  v6[0] = t6[0]; v6[1] = t6[1]; v6[2] = t6[2]; 
  v6[3] = t6[3]; v6[4] = t6[4]; v6[5] = t6[5]; 
}

void t6_to_v4

  ( const Vector& v4,
    const Vec6&   t6 )

{
  v4[0] = t6[0]; v4[1] = t6[1];
  v4[2] = t6[2]; v4[3] = t6[3];
}



//-----------------------------------------------------------------------
//   functions - tmatmul
//-----------------------------------------------------------------------

Vec6 tmatmul

  ( const Matrix& mat,
    const Vec6&   t6 )

{
  Vector v61_(6);
  Vector v62_(6);
  
  t6_to_v6 ( v61_, t6 );

  matmul ( v62_, mat, v61_ );

  return v6_to_t6 ( v62_ );
}


Vec6 tmatmul

  ( const Vec6&   t6,
    const Matrix& mat )

{
  Vector v61_(6);
  Vector v62_(6);
  
  t6_to_v6 ( v61_, t6 );

  matmul ( v62_, v61_, mat );

  return v6_to_t6 ( v62_ );

}

Vec6 tmatmul

  ( const Vec6&   t6,
    const Mat6&   tt6 )

{
  Matrix m6_1(6,6);
  tt6_to_m6 ( m6_1, tt6 );

  return tmatmul ( t6, m6_1 );
}

Vec6 tmatmul

  ( const Mat6&   tt6,
    const Vec6&   t6 )

{
  Matrix m6_1(6,6);
  tt6_to_m6 ( m6_1, tt6 );

  return tmatmul ( m6_1, t6 );
}


Mat6 tmatmul

  ( const Vec6&   t0,
    const Vec6&   t1  )

{
  Mat6   tt6;
  Vector v61_(6);
  Vector v62_(6);
  Matrix m6_ (6,6);
  
  t6_to_v6 ( v61_, t0 );
  t6_to_v6 ( v62_, t1 );
  
  matmul ( m6_, v61_, v62_ );

  m6_to_tt6 ( tt6,m6_ );

  return tt6;
}

Mat6 tmatmul

  ( const Mat6&     tt6,
    const Matrix&   mat  )
    
{
  Matrix m6_1(6,6);
  Matrix m6_2(6,6);
  Mat6   tmp6;
  
  tt6_to_m6 ( m6_1, tt6 );
  
  matmul ( m6_2 , m6_1, mat );

  m6_to_tt6 ( tmp6 , m6_2 );

  return tmp6;
}    
    
Mat6 tmatmul

  ( const Matrix&   mat,
    const Mat6&     tt6  )
        
{
  Matrix m6_1(6,6);
  Matrix m6_2(6,6);
  Mat6   tmp6;
  
  tt6_to_m6 ( m6_1, tt6 );
  
  matmul ( m6_2 , mat, m6_1 );

  m6_to_tt6 ( tmp6 , m6_2 );

  return tmp6;
}   

Mat6 tmatmul

  ( const Mat6&     tt6_1,
    const Mat6&     tt6_2  )
        
{
  Matrix m6_1(6,6);
  Matrix m6_2(6,6);
  Mat6   tmp6;
  
  tt6_to_m6 ( m6_1, tt6_1 );
  tt6_to_m6 ( m6_2, tt6_2 );

  m6_to_tt6 ( tmp6 , matmul ( m6_1, m6_2 ) );

  return tmp6;
}      
    
//-----------------------------------------------------------------------
//   functions - dot
//-----------------------------------------------------------------------

double tdot

  ( const Vector&   v6,
    const Vec6&     t6 )
        
{
  Vector v6_1(6);

  t6_to_v6( v6_1, t6 );
  
  return dot(v6,v6_1);
}   


double tdot

  ( const Vec6&     t6,
    const Vector&   v6 )
        
{
  Vector v6_1(6);

  t6_to_v6( v6_1, t6 );
  
  return dot(v6,v6_1);
}   