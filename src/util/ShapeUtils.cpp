/** 
 *  Some utility functions related to shape functions.
 *  Based on functions written by Vinh Phu Nguyen and
 *  Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
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
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/LUSolver.h>
#include <jem/util/Properties.h>
#include <jem/util/StringUtils.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/algebra/EigenUtils.h>
#include <jem/numeric/algebra/LUSolver.h>

#include <jive/geom/error.h>

#include "ShapeUtils.h"

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