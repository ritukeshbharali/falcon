
/** @file TensorUtils.cpp
 *  @brief Implements basic tensor operations.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */


#include <jem/base/Array.h>
#include <jive/util/utilities.h>
#include <jem/numeric/algebra/utilities.h>

#include "TensorUtils.h"

using namespace std;
using namespace jem;


//-----------------------------------------------------------------------------
//  trace
//-----------------------------------------------------------------------------

double tensorUtils::trace

  ( const Matrix&    a  )

{
  const int n = a.size(0);
  
  double c = 0.0;
  
  for (int i = 0; i < n; ++i)
  {
    c += a(i,i);
  }

  return c;  
}


//-----------------------------------------------------------------------------
//  otimes (c_ijkl = a_i b_j)
//-----------------------------------------------------------------------------

Matrix tensorUtils::otimes

  ( const Vector&    a ,
    const Vector&    b )

{
  const int arows    = a.size();
  const int brows    = b.size();
  
  Matrix c(arows, brows);
  c = 0.0;
  
  for (int i = 0; i < arows; ++i)
  {
    for (int j = 0; j < brows; ++j)
    {
      c(i,j) = a[i] * b[j];
    }
  }

  return c;  
}

//-----------------------------------------------------------------------------
//  otimes (c_ijkl = a_ij b_kl)
//-----------------------------------------------------------------------------

Quadix tensorUtils::otimes

  ( const Matrix&    a ,
    const Matrix&    b )

{
  const int arows    = a.size(0);
  const int acolumns = a.size(1);
  const int brows    = b.size(0);
  const int bcolumns = b.size(1);
  
  Quadix c(arows, acolumns, brows, bcolumns);
  c = 0.0;
  
  for (int i = 0; i < arows; ++i)
  {
    for (int j = 0; j < acolumns; ++j)
    {  
      for (int k = 0; k < brows; ++k)
      {
        for (int l = 0; l < bcolumns; ++l)
        {
          c(i,j,k,l) = a(i, j) * b(k, l);  
        }
      }
    }
  }

  return c;  
}

//------------------------------------------------------------------------------
//  otimesu (c_ijkl = a_ik b_jl)
//------------------------------------------------------------------------------

Quadix tensorUtils::otimesu

  ( const Matrix&    a ,
    const Matrix&    b )

{ 
  const int arows    = a.size(0);
  const int acolumns = a.size(1);
  const int brows    = b.size(0);
  const int bcolumns = b.size(1);
  
  Quadix c(arows, brows, acolumns, bcolumns);
  c = 0.0;
   
  for (int i = 0; i < arows; ++i)
  {
    for (int j = 0; j < brows; ++j)
    {
      for (int k = 0; k < acolumns; ++k)
      {
        for (int l = 0; l < bcolumns; ++l)
        {
          c(i,j,k,l) = a(i, k) * b(j, l);  
        }
      }
    }
  }

  return c;  
}

//------------------------------------------------------------------------------
//  otimesl (c_ijkl = a_il b_jk)
//------------------------------------------------------------------------------

Quadix tensorUtils::otimesl

  ( const Matrix&    a ,
    const Matrix&    b )

{ 
  const int arows    = a.size(0);
  const int acolumns = a.size(1);
  const int brows    = b.size(0);
  const int bcolumns = b.size(1);
  
  Quadix c(arows, brows, bcolumns, acolumns);
  c = 0.0;
  
  for (int i = 0; i < arows; ++i)
  {
    for (int j = 0; j < brows; ++j)
    {
      for (int k = 0; k < bcolumns; ++k)
      {
        for (int l = 0; l < acolumns; ++l)
        {
          c(i,j,k,l) = a(i, l) * b(j, k);  
        }
      }
    }
  } 
  return c;  
}

//------------------------------------------------------------------------------
//  doubleDot (c = a_ij b_ij)
//------------------------------------------------------------------------------

double tensorUtils::doubleDot

  ( const Matrix&    a ,
    const Matrix&    b )

{ 
  const int a_i      = a.size(0);
  const int a_j      = a.size(1);
  
  double c = 0.0;
  
  for (int i = 0; i < a_i; ++i)
  {  
    for (int j = 0; j < a_j; ++j)
    {
      c += a(i,j) * b(i,j);
    }
  }

  return c;  
}

//------------------------------------------------------------------------------
//  doubleDot (c_ij = a_ijkl b_kl) [overloaded]
//------------------------------------------------------------------------------

Matrix tensorUtils::doubleDot

  ( const Quadix&    a ,
    const Matrix&    b )

{ 
  const int a_i      = a.size(0);
  const int a_j      = a.size(1);
  const int brows    = b.size(0);
  const int bcolumns = b.size(1);
  
  Matrix c(a_i ,a_j );
  c = 0.0;
  
  for (int i = 0; i < a_i; ++i)
  {  
    for (int j = 0; j < a_j; ++j)
    {  
      for (int k = 0; k < brows; ++k)
      { 
        for (int l = 0; l < bcolumns; ++l)
        {  
          c(i,j) += a(i,j,k,l) * b(k,l);
        }
      }
    }
  }
  return c;  
}

// -----------------------------------------------------------------------
//   getI2 (second-order identity tensor)
// -----------------------------------------------------------------------

Matrix tensorUtils::getI2

  ( const int        rnk )

{

  Matrix A(rnk,rnk);
  A = 0.0;

  for ( int i = 0; i < rnk; i++ )
    for ( int j = 0; j < rnk; j++ )
      A(i,j) = (i==j);
  
  return A;
}

// -----------------------------------------------------------------------
//   getSymI4 (symmetric fourth-order identity tensor)
// -----------------------------------------------------------------------

Quadix tensorUtils::getSymI4

  ( const int        rnk )

{

  Quadix A(rnk,rnk,rnk,rnk);
  A = 0.0;

  for ( int i = 0; i < rnk; i++ )
    for ( int j = 0; j < rnk; j++ )
      for ( int k = 0; k < rnk; k++ )
        for ( int l = 0; l < rnk; l++ )
          A(i,j,k,l) = 0.5*((i==k)&&(j==l)) + 0.5*((i==l)&&(j==k));
  
  return A;
}

// -----------------------------------------------------------------------
//   getSymI4Dev (symmetric deviatoric fourth-order identity tensor)
// -----------------------------------------------------------------------

Quadix tensorUtils::getSymI4Dev

  ( const int        rnk )

{

  Quadix A(rnk,rnk,rnk,rnk);
  A = 0.0;

  for ( int i = 0; i < rnk; i++ )
    for ( int j = 0; j < rnk; j++ )
      for ( int k = 0; k < rnk; k++ )
        for ( int l = 0; l < rnk; l++ )
        {
          A(i,j,k,l) =     0.5*((i==k)&&(j==l)) + 0.5*((i==l)&&(j==k)) 
                       - 1./3.*((i==j)&&(k==l)) ;
        }
  
  return A;
}

// -----------------------------------------------------------------------
//   tensor2voigt (symmetrized)
// -----------------------------------------------------------------------

Matrix tensorUtils::tensor2voigt

    (  const Quadix&   A,
       const int       strCount )

{

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Matrix M(strCount,strCount);
  M = 0.0;

  if ( strCount == 4 )   // 2D {xx, yy, zz, xy}
  {
    M(0,0) = A(0,0,0,0);                          // xx-xx
    M(1,0) = A(1,1,0,0);                          // yy-xx
    M(2,0) = A(2,2,0,0);                          // zz-xx
    M(3,0) = 0.5 * ( A(0,1,0,0) + A(1,0,0,0) );   // xy-xx

    M(0,1) = A(0,0,1,1);                          // xx-yy
    M(1,1) = A(1,1,1,1);                          // yy-yy
    M(2,1) = A(2,2,1,1);                          // zz-yy
    M(3,1) = 0.5 * ( A(0,1,1,1) + A(1,0,1,1) );   // xy-yy

    M(0,2) = A(0,0,2,2);                          // xx-zz
    M(1,2) = A(1,1,2,2);                          // yy-zz
    M(2,2) = A(2,2,2,2);                          // zz-zz
    M(3,2) = 0.5 * ( A(0,1,2,2) + A(1,0,2,2) );   // xy-zz

    M(0,3) = 0.5 * ( A(0,0,0,1) + A(0,0,1,0) );   // xx-xy
    M(1,3) = 0.5 * ( A(1,1,0,1) + A(1,1,1,0) );   // yy-xy
    M(2,3) = 0.5 * ( A(2,2,0,1) + A(2,2,1,0) );   // zz-xy
    M(3,3) = 0.25 *( A(0,1,0,1) + A(1,0,0,1) 
                   + A(0,1,1,0) + A(1,0,1,0) );   // xy-xy

  }

  else if ( strCount == 6 )   // 3D {xx, yy, zz, xy, yz, zx}
  {

    M(0,0) = A(0,0,0,0);                          // xx-xx
    M(1,0) = A(1,1,0,0);                          // yy-xx
    M(2,0) = A(2,2,0,0);                          // zz-xx
    M(3,0) = 0.5 * ( A(0,1,0,0) + A(1,0,0,0) );   // xy-xx
    M(4,0) = 0.5 * ( A(1,2,0,0) + A(2,1,0,0) );   // yz-xx
    M(5,0) = 0.5 * ( A(2,0,0,0) + A(0,2,0,0) );   // zx-xx

    M(0,1) = A(0,0,1,1);                          // xx-yy
    M(1,1) = A(1,1,1,1);                          // yy-yy
    M(2,1) = A(2,2,1,1);                          // zz-yy
    M(3,1) = 0.5 * ( A(0,1,1,1) + A(1,0,1,1) );   // xy-yy
    M(4,1) = 0.5 * ( A(1,2,1,1) + A(2,1,1,1) );   // yz-yy
    M(5,1) = 0.5 * ( A(2,0,1,1) + A(0,2,1,1) );   // zx-yy

    M(0,2) = A(0,0,2,2);                          // xx-zz
    M(1,2) = A(1,1,2,2);                          // yy-zz
    M(2,2) = A(2,2,2,2);                          // zz-zz
    M(3,2) = 0.5 * ( A(0,1,2,2) + A(1,0,2,2) );   // xy-zz
    M(4,2) = 0.5 * ( A(1,2,2,2) + A(2,1,2,2) );   // yz-zz
    M(5,2) = 0.5 * ( A(2,0,2,2) + A(0,2,2,2) );   // zx-zz

    M(0,3) = 0.5 * ( A(0,0,0,1) + A(0,0,1,0) );   // xx-xy
    M(1,3) = 0.5 * ( A(1,1,0,1) + A(1,1,1,0) );   // yy-xy
    M(2,3) = 0.5 * ( A(2,2,0,1) + A(2,2,1,0) );   // zz-xy
    M(3,3) = 0.25 *( A(0,1,0,1) + A(1,0,1,0) 
                   + A(0,1,0,1) + A(1,0,1,0) );   // xy-xy
    M(4,3) = 0.25 *( A(1,2,0,1) + A(2,1,1,0) 
                   + A(2,1,0,1) + A(1,2,1,0) );   // yz-xy
    M(5,3) = 0.25 *( A(0,2,0,1) + A(2,0,1,0) 
                   + A(2,0,0,1) + A(0,2,1,0) );   // zx-xy

    M(0,4) = 0.5 * ( A(0,0,1,2) + A(0,0,2,1) );   // xx-yz
    M(1,4) = 0.5 * ( A(1,1,1,2) + A(1,1,2,1) );   // yy-yz
    M(2,4) = 0.5 * ( A(2,2,1,2) + A(2,2,2,1) );   // zz-yz
    M(3,4) = 0.25 *( A(0,1,1,2) + A(1,0,2,1) 
                   + A(0,1,1,2) + A(1,0,2,1) );   // xy-yz
    M(4,4) = 0.25 *( A(1,2,1,2) + A(2,1,2,1) 
                   + A(2,1,1,2) + A(1,2,2,1) );   // yz-yz
    M(5,4) = 0.25 *( A(0,2,1,2) + A(2,0,2,1) 
                   + A(2,0,1,2) + A(0,2,2,1) );   // zx-yz

    M(0,5) = 0.5 * ( A(0,0,0,2) + A(0,0,2,0) );   // xx-zx
    M(1,5) = 0.5 * ( A(1,1,0,2) + A(1,1,2,0) );   // yy-zx
    M(2,5) = 0.5 * ( A(2,2,0,2) + A(2,2,2,0) );   // zz-zx
    M(3,5) = 0.25 *( A(0,1,0,2) + A(1,0,2,0) 
                   + A(0,1,0,2) + A(1,0,2,0) );   // xy-zx
    M(4,5) = 0.25 *( A(1,2,0,2) + A(2,1,2,0) 
                   + A(2,1,0,2) + A(1,2,2,0) );   // yz-zx
    M(5,5) = 0.25 *( A(0,2,0,2) + A(2,0,2,0) 
                   + A(2,0,0,2) + A(0,2,2,0) );   // zx-zx

  }
  
  return M;
}


// -----------------------------------------------------------------------
//   tensor2voigtStrain
// -----------------------------------------------------------------------

Vector tensorUtils::tensor2voigtStrain

    (  const Matrix&   A,
       const int       strCount )
{

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Vector c (strCount);
  c = 0.0;

  c[0] = A(0,0);
  c[1] = A(1,1);
  c[2] = A(2,2);
  c[3] = A(0,1) + A(1,0);

  if ( strCount == 6 )
  {
    c[4] = A(1,2) + A(2,1);
    c[5] = A(0,2) + A(2,0);
  }

  return c;

}


// -----------------------------------------------------------------------
//   tensor2voigtStress
// -----------------------------------------------------------------------

Vector tensorUtils::tensor2voigtStress

    (  const Matrix&   A,
       const int       strCount )
{

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Vector c (strCount);
  c = 0.0;

  c[0] = A(0,0);
  c[1] = A(1,1);
  c[2] = A(2,2);
  c[3] = 0.5 * ( A(0,1) + A(1,0) );

  if ( strCount == 6 )
  {
    c[4] = 0.5 * ( A(1,2) + A(2,1) );
    c[5] = 0.5 * ( A(0,2) + A(2,0) );
  }

  return c;
  
}


// -----------------------------------------------------------------------
//   voigt2tensorStrain
// -----------------------------------------------------------------------

Matrix tensorUtils::voigt2tensorStrain

    ( const Vector&   c )
{

  const int strCount = c.size();

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Matrix A (3,3);
  A = 0.0;

  A(0,0) = c[0];
  A(1,1) = c[1];
  A(2,2) = c[2];
  A(0,1) = A(1,0) = .5 * c[3];

  if ( strCount == 6 )
  {
    A(1,2) = A(2,1) = .5 * c[4];
    A(2,0) = A(0,2) = .5 * c[5];
  }

  return A;
  
}


// -----------------------------------------------------------------------
//   voigt2tensorStress
// -----------------------------------------------------------------------

Matrix tensorUtils::voigt2tensorStress

    ( const Vector&   c )
{

  const int strCount = c.size();

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Matrix A (3,3);
  A = 0.0;

  A(0,0) = c[0];
  A(1,1) = c[1];
  A(2,2) = c[2];
  A(0,1) = A(1,0) = c[3];

  if ( strCount == 6 )
  {
    A(1,2) = A(2,1) = c[4];
    A(2,0) = A(0,2) = c[5];
  }

  return A;
  
}


// -----------------------------------------------------------------------
//   voigt2tensorRankStrain
// -----------------------------------------------------------------------

Matrix tensorUtils::voigt2tensorRankStrain

    ( const Vector&   c )
{

  const int strCount = c.size();

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Matrix A;

  if ( strCount == 4 )
  {
    A.resize (2,2);
    A = 0.0;

   A(0,0) = c[0];
   A(1,1) = c[1];
   A(0,1) = A(1,0) = .5 * c[3];
  }
  else
  {
    A.resize (3,3);
    A = 0.0;

   A(0,0) = c[0];
   A(1,1) = c[1];
   A(2,2) = c[2];
   A(0,1) = A(1,0) = .5 * c[3];
   A(1,2) = A(2,1) = .5 * c[4];
   A(2,0) = A(0,2) = .5 * c[5];

  } 

  return A;
  
}


// -----------------------------------------------------------------------
//   voigt2tensorRankStress
// -----------------------------------------------------------------------

Matrix tensorUtils::voigt2tensorRankStress

    ( const Vector&   c )
{

  const int strCount = c.size();

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Matrix A;

  if ( strCount == 4 )
  {
    A.resize (2,2);
    A = 0.0;

   A(0,0) = c[0];
   A(1,1) = c[1];
   A(0,1) = A(1,0) = c[3];
  }
  else
  {
    A.resize (3,3);
    A = 0.0;

   A(0,0) = c[0];
   A(1,1) = c[1];
   A(2,2) = c[2];
   A(0,1) = A(1,0) = c[3];
   A(1,2) = A(2,1) = c[4];
   A(2,0) = A(0,2) = c[5];

  } 

  return A;
  
}


// -----------------------------------------------------------------------
//   engng2tensorial
// -----------------------------------------------------------------------

Vector tensorUtils::engng2tensorial

    ( const Vector&   c )
{

  const int strCount = c.size();

  JEM_ASSERT ( strCount == 4 || strCount == 6 );

  Vector d (strCount);
  d = c;

  d[3] *= 0.5;

  if ( strCount == 6 )
  {
    d[4] *= 0.5;
    d[5] *= 0.5;
  }

  return d;
  
}


// ---------------------------------------------------------------------
//  fill2Dstress (3D to 2D)
// ---------------------------------------------------------------------

Vector tensorUtils::fill2Dstress

   ( const Vector&    a  )
  
{
  Vector b (4);

  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  b[3] = a[3];

  return b;
}


// ---------------------------------------------------------------------
//  fill2Dstrain (3D to 2D)
// ---------------------------------------------------------------------

Vector tensorUtils::fill2Dstrain

   ( const Vector&    a  )
  
{
  Vector b (4);

  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  b[3] = a[3];

  return b;
}


// ---------------------------------------------------------------------
//  fill3Dstress (2D to 3D)
// ---------------------------------------------------------------------

Vector tensorUtils::fill3Dstress 

   ( const Vector&    a  )
  
{
  Vector b (6);

  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  b[3] = a[3];
  b[4] = 0.0;
  b[5] = 0.0;    

  return b;
}


// ---------------------------------------------------------------------
//  fill3Dtensorialstrain (2D to 3D)
// ---------------------------------------------------------------------

Vector tensorUtils::fill3Dtensorialstrain

   ( const Vector&    a  )
  
{
  Vector b (6);

  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  b[3] = 0.5 * a[3];
  b[4] = 0.0;
  b[5] = 0.0;    

  return b;
}


// ---------------------------------------------------------------------
//  fill3Dstrain (2D to 3D)
// ---------------------------------------------------------------------

Vector tensorUtils::fill3Dstrain

   ( const Vector&    a  )
  
{
  Vector b (6);

  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
  b[3] = a[3];
  b[4] = 0.0;
  b[5] = 0.0;    

  return b;
}


// ---------------------------------------------------------------------
//  fill2Dmatrix (3D to 2D)
// ---------------------------------------------------------------------

Matrix tensorUtils::fill2Dmatrix

   ( const Matrix&    A  )
  
{
  Matrix B (4,4);

  B(0,0) = A(0,0);
  B(0,1) = A(0,1);
  B(0,2) = A(0,2);
  B(0,3) = A(0,3);

  B(1,0) = A(1,0);
  B(1,1) = A(1,1);
  B(1,2) = A(1,2);
  B(1,3) = A(1,3);

  B(2,0) = A(2,0);
  B(2,1) = A(2,1);
  B(2,2) = A(2,2);
  B(2,3) = A(2,3);

  B(3,0) = A(3,0);
  B(3,1) = A(3,1);
  B(3,2) = A(3,2);
  B(3,3) = A(3,3);

  return B;
}


