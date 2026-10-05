/** 
 *  @file ShapeUtils.cpp
 *  @brief Implementation of utility functions related to shape functions, 
 *  strain-displacement matrices, and finite strain kinematics.
 *  Based on functions written by Vinh Phu Nguyen and Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [05 October 2026] Removed functions not needed. Added functions
 *       for large deformation/finite strain. (RB)
 */

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>

#include <jem/io/PrintWriter.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/System.h>
#include <jem/base/Error.h>
#include <jem/base/Float.h>
#include <jem/base/array/tensor.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/EigenUtils.h>
#include <jem/numeric/algebra/LUSolver.h>

#include <jive/geom/error.h>

#include "ShapeUtils.h"
#include "Constants.h"
#include "TensorUtils.h"

using jem::System;
using jem::ALL;
using jem::END;
using jem::TensorIndex;
using jem::io::endl;
using jem::idx_t;
using jem::Array;

//=======================================================================
//   B matrix
//=======================================================================

//-----------------------------------------------------------------------
//   get1DShapeGrads
//-----------------------------------------------------------------------


void              get1DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_PRECHECK ( b.size(0) == 1 &&
                 g.size(0) == 1 &&
                 b.size(1) == g.size(1) );

  b = g;
}


//-----------------------------------------------------------------------
//   get2DShapeGrads
//-----------------------------------------------------------------------

// VP Nguyen, 6 October 2014
// in 2D, B matrix has dimension 4x2n where n is the number of nodes
// epsilon_xy in the last row. This is needed for constitutive models
// where sigma_zz and epsilon_zz are required for 2D problems.

void              get2DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_PRECHECK ( b.size(0) == 4 &&
                 g.size(0) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  int    i, i1;
  double Nix, Niy;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    i  = 2 * inode;
    i1 = i + 1;

    Nix = g(0,inode);
    Niy = g(1,inode);

    b(0,i ) = Nix;
    b(1,i1) = Niy;

    b(3,i ) = Niy;
    b(3,i1) = Nix;
  }
}


//-----------------------------------------------------------------------
//   get3DShapeGrads
//-----------------------------------------------------------------------

// The strain-displacement matrix B is for a strain vector stored
// as [epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy, epsilon_yz, epsilon_zx].

void              get3DShapeGrads

  ( const Matrix&   b,
    const Matrix&   g )

{
  JEM_PRECHECK ( b.size(0) == 6 &&
                 g.size(0) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 3 * inode;

    b(0,i + 0) = g(0,inode);
    b(1,i + 1) = g(1,inode);
    b(2,i + 2) = g(2,inode);

    b(3,i + 0) = g(1,inode);
    b(3,i + 1) = g(0,inode);

    b(4,i + 1) = g(2,inode);
    b(4,i + 2) = g(1,inode);

    b(5,i + 2) = g(0,inode);
    b(5,i + 0) = g(2,inode);
  }
}


//-----------------------------------------------------------------------
//   getShapeGradsFunc
//-----------------------------------------------------------------------


ShapeGradsFunc getShapeGradsFunc ( int rank )
{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );

  if      ( rank == 1 )
  {
    return & get1DShapeGrads;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeGrads;
  }
  else
  {
    return & get3DShapeGrads;
  }
}


//=======================================================================
//   B matrix (Linear): used for Finite Strain problems
//=======================================================================

//-----------------------------------------------------------------------
//    getBMatrixLin2D
//-----------------------------------------------------------------------

void  getBMatrixLin2D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g )

{
  JEM_ASSERT   ( b.size(0) == 4 &&
                 g.size(0) == 2 &&
                 f.size(0) == 2 &&
                 f.size(1) == 2 &&
                 b.size(1) == 2 * g.size(1) );

  const int  nodeCount = g.size (1);

  b = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i0 = 2 * inode;
    int  i1 = i0 + 1;

    b(0,i0) = f(0,0) * g(0,inode);
    b(0,i1) = f(1,0) * g(0,inode);

    b(1,i0) = f(0,1) * g(1,inode);
    b(1,i1) = f(1,1) * g(1,inode);

    b(3,i0) = f(0,0) * g(1,inode) + f(0,1) * g(0,inode);
    b(3,i1) = f(1,1) * g(0,inode) + f(1,0) * g(1,inode);
  }
}

//-----------------------------------------------------------------------
//    getBMatrixLin3D
//-----------------------------------------------------------------------

void  getBMatrixLin3D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g )

{
  JEM_ASSERT   ( b.size(0) == 6 &&
                 g.size(0) == 3 &&
                 f.size(0) == 3 &&
                 f.size(1) == 3 &&
                 b.size(1) == 3 * g.size(1) );

  const int  nodeCount = g.size (1);

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i0 = 3 * inode;
    int  i1 = i0 + 1;
    int  i2 = i0 + 2;

    b(0,i0) = f(0,0) * g(0,inode);
    b(0,i1) = f(1,0) * g(0,inode);
    b(0,i2) = f(2,0) * g(0,inode);

    b(1,i0) = f(0,1) * g(1,inode);
    b(1,i1) = f(1,1) * g(1,inode);
    b(1,i2) = f(2,1) * g(1,inode);

    b(2,i0) = f(0,2) * g(2,inode);
    b(2,i1) = f(1,2) * g(2,inode);
    b(2,i2) = f(2,2) * g(2,inode);

    b(3,i0) = f(0,0) * g(1,inode) + f(0,1) * g(0,inode);
    b(3,i1) = f(1,0) * g(1,inode) + f(1,1) * g(0,inode);
    b(3,i2) = f(2,0) * g(1,inode) + f(2,1) * g(0,inode);

    b(4,i0) = f(0,1) * g(2,inode) + f(0,2) * g(1,inode);
    b(4,i1) = f(1,1) * g(2,inode) + f(1,2) * g(1,inode);
    b(4,i2) = f(2,1) * g(2,inode) + f(2,2) * g(1,inode);

    b(5,i0) = f(0,2) * g(0,inode) + f(0,0) * g(2,inode);
    b(5,i1) = f(1,2) * g(0,inode) + f(1,0) * g(2,inode);
    b(5,i2) = f(2,2) * g(0,inode) + f(2,0) * g(2,inode);
  }
}

//-----------------------------------------------------------------------
//   getBMatrixLinFunc
//-----------------------------------------------------------------------

BMatrixLinFunc getBMatrixLinFunc ( int rank )
{
  JEM_PRECHECK ( rank > 1 && rank <= 3 );

  if ( rank == 2 )
  {
    return & getBMatrixLin2D;
  }
  else
  {
    return & getBMatrixLin3D;
  }
}


//=======================================================================
//   N matrix
//=======================================================================

// --------------------------------------------------------------------
//  get1DShapeFuncs
// --------------------------------------------------------------------

void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n )
{
  sfuncs(0,0) = n[0];
}

// --------------------------------------------------------------------
//  get2DShapeFuncs
// --------------------------------------------------------------------

void                  get2DShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_PRECHECK ( s.size(0) == 2 &&
                 s.size(1) == 2 * n.size() );

  const int  nodeCount = n.size ();

  s = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 2 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
  }
}

// --------------------------------------------------------------------
//  get3DShapeFuncs
// --------------------------------------------------------------------

void                  get3DShapeFuncs

  ( const Matrix&       s,
    const Vector&       n )
{
  JEM_PRECHECK ( s.size(0) == 3 &&
                 s.size(1) == 3 * n.size() );

  const int  nodeCount = n.size ();

  s = 0.0;

  for ( int inode = 0; inode < nodeCount; inode++ )
  {
    int  i = 3 * inode;

    s(0,i + 0) = n[inode];
    s(1,i + 1) = n[inode];
    s(2,i + 2) = n[inode];
  }
}

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeFunc        getShapeFunc

  ( int                 rank )

{
  JEM_PRECHECK ( rank >= 1 && rank <= 3 );


  if      ( rank == 1 )
  {
    return & get1DShapeFuncs;
  }
  else if ( rank == 2 )
  {
    return & get2DShapeFuncs;
  }
  else
  {
    return & get3DShapeFuncs;
  }
}

//=======================================================================
//   Finite Strain utility functions
//=======================================================================

//-----------------------------------------------------------------------
//  Add geometric nonlinear part to the stiffness matrix k. This requires
//  the Second Piola Kirchoff stress pk2 and the shape function gradients
//  g and the weight function w.
//-----------------------------------------------------------------------

void  addGeomNonlinElemMat

  (       Matrix&  k,
    const Vector&  pk2,
    const Matrix&  g,
    const double   w )
{
  int rank = g.size(0);
  int nn   = g.size(1);

  JEM_ASSERT    ( k.size(0) == rank*nn &&
                  k.size(1) == rank*nn );

  Matrix   t = tensorUtils::voigt2tensorRankStress( pk2 );

  JEM_ASSERT    ( t.size(0) == rank &&
                  t.size(1) == rank );

  Matrix   btb ( nn, nn );

  TensorIndex in,jn,it,jt;

  btb(in,jn) = dot( g(jt,in), dot( t(it,jt), g(it,jn), it ), jt );

  for ( int in = 0; in < nn; ++in )
  {
    for ( int jn = 0; jn < nn; ++jn )
    {
      for ( int ix = 0; ix < rank; ++ix )
      {
        k( in*rank+ix, jn*rank+ix ) += w*btb(in,jn);
      }
    }
  }
}

void  addGeomNonlinElemMat

  (       Matrix&  k,
    const Matrix&  pk2,
    const Matrix&  g,
    const double   w )
{
  int rank = g.size(0);
  int nn   = g.size(1);

  JEM_ASSERT    ( k.size(0) == rank*nn &&
                  k.size(1) == rank*nn );

  JEM_ASSERT    ( pk2.size(0) == rank &&
                  pk2.size(1) == rank );

  Matrix   btb ( nn, nn );

  TensorIndex in,jn,it,jt;

  btb(in,jn) = dot( g(jt,in), dot( pk2(it,jt), g(it,jn), it ), jt );

  for ( int in = 0; in < nn; ++in )
  {
    for ( int jn = 0; jn < nn; ++jn )
    {
      for ( int ix = 0; ix < rank; ++ix )
      {
        k( in*rank+ix, jn*rank+ix ) += w*btb(in,jn);
      }
    }
  }

}

// -----------------------------------------------------------------------
//   Deformation gradient
// -----------------------------------------------------------------------

void                  getDeformationGradient

  (       Matrix&       f,
    const Vector&       u,
    const Matrix&       g )
{
  idx_t rank = g.size(0);

  for (idx_t i = 0; i < rank; ++i)
  {
    for (idx_t j = 0; j < rank; ++j)
    {
      f(i,j) = dot ( g(j,ALL), u[slice(i,END,rank)] );
    }
    f(i,i) += 1.0;
  }
}


// -----------------------------------------------------------------------
//   Green-Lagrange strain
// -----------------------------------------------------------------------

void                  getGreenLagrangeStrain

  (       Matrix&       e,
    const Matrix&       f )
{
  idx_t rank = f.size(0);

  // Compute Green-Lagrange tensor e = 0.5 (F^T F - I)

  TensorIndex i,j,k;
  e.resize (rank, rank);

  e(i,j) = 0.5 * ( dot( f(k,i), f(k,j), k ) - where(i==j,1.0,0.0) );

}

// -----------------------------------------------------------------------
//   Green-Lagrange strain
// -----------------------------------------------------------------------

void                  getGreenLagrangeStrain

  (       Vector&       eps,
    const Matrix&       f )
{
  idx_t rank = f.size(0);

  Matrix e (rank, rank);

  getGreenLagrangeStrain( e, f);

  // Convert to Voigt notation

  eps = tensorUtils::tensor2voigtStrain( e, STRAIN_COUNTS[rank] );
}

// -----------------------------------------------------------------------
//   1st Piola stress
// -----------------------------------------------------------------------

void                  getFirstPiolaStress

  (       Matrix&       pk1,
    const Matrix&       f,
    const Matrix&       pk2 )
{
  int rank = f.size(0);

  // Ensure strict dimension compatibility for (rank x rank) matrices:
  // - F must be square (rank x rank)
  // - PK2 must be square (rank x rank)
  // - PK1 must match (rank x rank)
  JEM_ASSERT( f.size(1)   == rank &&
              pk2.size(0) == rank &&
              pk2.size(1) == rank &&
              pk1.size(0) == rank &&
              pk1.size(1) == rank );

  MChain2 mc2;

  pk1 = mc2.matmul( f, pk2 );
}


void                  getFirstPiolaStress

  (       Matrix&       pk1,
    const Matrix&       f,
    const Vector&       pk2 )
{
  int rank = f.size(0);

  // Ensure strict dimension compatibility for (rank x rank) matrices:
  // - F must be square (rank x rank)
  JEM_ASSERT( f.size(1)   == rank );

  Matrix pk2_mat = tensorUtils::voigt2tensorRankStress( pk2 );

  JEM_ASSERT( pk2_mat.size(0)   == rank &&
              pk2_mat.size(1)   == rank );

  MChain2 mc2;

  pk1 = mc2.matmul( f, pk2_mat );
}

// -----------------------------------------------------------------------
//   Cauchy stress
// -----------------------------------------------------------------------

void                  getCauchyStress

  (       Vector&       sigma,
    const Matrix&       f,
    const Matrix&       pk2 )
{
  idx_t rank = f.size(0);

  JEM_ASSERT ( f.size(1)   == rank &&
               pk2.size(0) == rank &&
               pk2.size(1) == rank );

  MChain3 mc3;

  const double J = jem::numeric::det(f);

  Matrix sigma_mat (rank, rank);

  sigma_mat = (mc3.matmul( f, pk2, f.transpose() ) ) / J;

  sigma = tensorUtils::tensor2voigtStress( sigma_mat,
                                           STRAIN_COUNTS[rank] );
  
}


void getCauchyStress
(
        Vector& sigma,
  const Matrix& f,
  const Vector& pk2
)
{
  const idx_t rank = f.size(0);

  MChain3 mc3;

  if ( rank == 2 )
  {
    // Embed 2D deformation gradient in 3D
    Matrix f3(3,3);
    f3 = 0.0;

    f3(0,0) = f(0,0);
    f3(0,1) = f(0,1);
    f3(1,0) = f(1,0);
    f3(1,1) = f(1,1);
    f3(2,2) = 1.0;

    // Full {xx,yy,zz,xy} PK2 -> 3x3 tensor
    Matrix pk2_mat =
        tensorUtils::voigt2tensorStress(pk2);

    const double J = jem::numeric::det(f);

    Matrix sigma_mat(3,3);

    sigma_mat =
        mc3.matmul(f3, pk2_mat, f3.transpose()) / J;

    sigma =
        tensorUtils::tensor2voigtStress(sigma_mat,4);
  }
  else
  {
    Matrix pk2_mat =
        tensorUtils::voigt2tensorStress(pk2);

    const double J = jem::numeric::det(f);

    Matrix sigma_mat(3,3);

    sigma_mat =
        mc3.matmul(f, pk2_mat, f.transpose()) / J;

    sigma =
        tensorUtils::tensor2voigtStress(sigma_mat,6);
  }
}