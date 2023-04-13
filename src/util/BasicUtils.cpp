/** 
 *  Some basic functions common to most FE models and
 *  post-processing operations are implemented. Almost
 *  all functions are built on top of codes by written
 *  by Vinh Phu Nguyen and Frans van der Meer.
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

#include "BasicUtils.h"

#define tolerance 0.1e-20

extern "C"
{
  #include  <math.h>
}

using jem::System;
using jem::ALL;
using jem::END;
using jem::TensorIndex;
using jem::io::endl;
using jem::idx_t;
using jem::Array;


//-----------------------------------------------------------------------
//   constants
//-----------------------------------------------------------------------

// number of strain components (Voigt notation)
// corresponding to the rank of the problem
//                                  1D  2D  3D
const int    STRAIN_COUNTS[4] = { 0, 1, 4, 6 };

const double PI               = 3.14159265;
const double PI3              = 0.3333333333 * PI;
const double ONE_THIRD        = 0.33333333333;
const double TWO_THIRD        = 0.66666666667;
const double EPS              = 1.e-16;
const double RAD_120          = 0.6666666666666667 * PI;
const double RAD_240          = 2.0 * RAD_120;


//-----------------------------------------------------------------------
//   Bubble shape functions and derivatives (Triangles)
//-----------------------------------------------------------------------

Matrix                getBubbleShapeFunctions

  ( const Matrix&       ischeme,
    const Matrix&       coords  )

{
  // Only Triangle3 with minimum 3 int points
  JEM_PRECHECK ( ischeme.size(0) == 3 &&
                 ischeme.size(1) > 2  );

  Matrix sfuncs ( 4, ischeme.size(1) );

  sfuncs = 0.0;

  for ( int i = 0; i < ischeme.size(1); i++ )
  {
    const double xi  = ischeme(1,i);
    const double eta = ischeme(2,i);

    sfuncs( 3, i )   = xi * eta * ( 1. - xi - eta ) ;

    sfuncs( 0, i )   = ( 1. - eta - xi ) - 0.3333 * sfuncs( 3, i );
    sfuncs( 1, i )   = eta               - 0.3333 * sfuncs( 3, i );
    sfuncs( 2, i )   = xi                - 0.3333 * sfuncs( 3, i );

  }

  return sfuncs;
}  


Cubix               getBubbleShapeGradients

  ( const Matrix&       ischeme,
    const Matrix&       coords  )

{
  // Only Triangle3 with minimum 3 int points
  JEM_PRECHECK ( ischeme.size(0) == 3 &&
                 ischeme.size(1) > 2  );

  //           rank, nodeCount, ipCount              
  Cubix  grads  ( 2, 4, ischeme.size(1) );
  Matrix invJac ( 2, 2);

  grads  = 0.0;
  invJac = 0.0;

  // Fill Jac

  invJac(0,0) = coords(1,0) - coords(0,0);
  invJac(0,1) = coords(1,1) - coords(0,1);
  invJac(1,0) = coords(2,0) - coords(0,0);
  invJac(1,1) = coords(2,1) - coords(0,1);

  // Invert Jac

  using jem::numeric::LUSolver;

  double d;

  LUSolver::invert ( invJac, d );

  for ( int i = 0; i < ischeme.size(1); i++ )
  {
    const double xi  = ischeme(1,i);
    const double eta = ischeme(2,i);

    const double dN4dxi  = (eta * ( 1. - xi - eta ) - xi * eta);
    const double dN4deta = (xi  * ( 1. - xi - eta ) - xi * eta);

    const double dN1dxi  = -1. - 0.3333 * dN4dxi;
    const double dN1deta = -1. - 0.3333 * dN4deta;

    const double dN2dxi  = 1. - 0.3333 * dN4dxi;
    const double dN2deta = 0. - 0.3333 * dN4deta;

    const double dN3dxi  = 0. - 0.3333 * dN4dxi;
    const double dN3deta = 1. - 0.3333 * dN4deta;

    grads( 0, 3, i )    =   dN4dxi  * invJac(0,0) 
                          + dN4deta * invJac(0,1);

    grads( 1, 3, i )    =   dN4dxi  * invJac(1,0) 
                          + dN4deta * invJac(1,1);

    grads( 0, 0, i )    =   dN1dxi  * invJac(0,0) 
                          + dN1deta * invJac(0,1);

    grads( 1, 0, i )    =   dN1dxi  * invJac(1,0) 
                          + dN1deta * invJac(1,1);

    grads( 0, 1, i )    =   dN2dxi  * invJac(0,0) 
                          + dN2deta * invJac(0,1);

    grads( 1, 1, i )    =   dN2dxi  * invJac(1,0) 
                          + dN2deta * invJac(1,1);                      

    grads( 0, 2, i )    =   dN3dxi  * invJac(0,0) 
                          + dN3deta * invJac(0,1);

    grads( 1, 2, i )    =   dN3dxi  * invJac(1,0) 
                          + dN3deta * invJac(1,1);
    
  }

  return grads;
}


//-----------------------------------------------------------------------
//   Bubble shape functions and derivatives (Quads)
//-----------------------------------------------------------------------

Matrix                getQuadBubbleShapeFunctions

  ( const Matrix&       ischeme,
    const Matrix&       coords  )

{
  // Only Triangle3 with 6 int points
  /*JEM_PRECHECK ( ischeme.size(0) == 3 &&
                 ischeme.size(1) > 5  );*/

  Matrix sfuncs ( 5, ischeme.size(1) );

  sfuncs = 0.0;

  for ( int i = 0; i < ischeme.size(1); i++ )
  {
    const double xi  = ischeme(1,i);
    const double eta = ischeme(2,i);

    sfuncs( 4, i )   = ( 1. - xi - xi ) * ( 1. - eta - eta ) ;

  }

  return sfuncs;
}  


Cubix               getQuadBubbleShapeGradients

  ( const Matrix&       ischeme,
    const Matrix&       coords  )

{
  // Only Triangle3 with 6 int points
  /*JEM_PRECHECK ( ischeme.size(0) == 3 &&
                 ischeme.size(1) > 5  );*/

  //           rank, nodeCount, ipCount              
  Cubix  grads  ( 2, 5, ischeme.size(1) );
  Matrix invJac ( 2, 2);

  grads  = 0.0;
  invJac = 0.0;

  // Fill Jac

  invJac(0,0) = coords(1,0) - coords(0,0);
  invJac(0,1) = coords(1,1) - coords(0,1);
  invJac(1,0) = coords(2,0) - coords(0,0);
  invJac(1,1) = coords(2,1) - coords(0,1);

  // Invert Jac

  using jem::numeric::LUSolver;

  double d;

  LUSolver::invert ( invJac, d );

  for ( int i = 0; i < ischeme.size(1); i++ )
  {
    const double xi  = ischeme(1,i);
    const double eta = ischeme(2,i);

    const double dN5dxi  = - 2. * xi  * ( 1. - eta - eta );
    const double dN5deta = - 2. * eta * ( 1. - xi  - xi  );

    grads( 0, 4, i )    =   dN5dxi  * invJac(0,0) 
                          + dN5deta * invJac(0,1);

    grads( 1, 4, i )    =   dN5dxi  * invJac(1,0) 
                          + dN5deta * invJac(1,1);

  }

  return grads;
}


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


// --------------------------------------------------------------------
//  get2DNormalFuncs
// --------------------------------------------------------------------

void                  get2DNormalFuncs

  ( const Matrix&       normal,
    const Vector&       normalV )
{
   normal(0,0) = normalV[0]; normal(0,2) = normalV[1];
   normal(1,1) = normalV[1]; normal(1,2) = normalV[0];

}

// --------------------------------------------------------------------
//  get3DNormalFuncs
// --------------------------------------------------------------------

void                  get3DNormalFuncs

  ( const Matrix&       normal,
    const Vector&       normalV )
{
   normal(0,0) = normalV[0]; normal(0,3) = normalV[1]; normal(0,5) = normalV[2];
   normal(1,1) = normalV[1]; normal(1,3) = normalV[0]; normal(1,4) = normalV[2];
   normal(2,2) = normalV[2]; normal(2,4) = normalV[1]; normal(2,5) = normalV[0];

}

// A function that returns a pointer to a function that computes the
// matrix of normal functions given the number of spatial dimensions.

NormalMatrixFunc        getNormalFunc

  ( int                 rank )

{
  if ( rank == 2 )
  {
    return & get2DNormalFuncs;
  }
  else
  {
    return & get3DNormalFuncs;
  }
}

void                  get2DNormalInMatrixFuncs

  ( const Matrix&       normal,
    const Vector&       normalV )
{
   normal(0,0) = normalV[0]; normal(0,3) = normalV[1];
   normal(1,1) = normalV[1]; normal(1,3) = normalV[0];

}

// A function that returns a pointer to a function that computes the
// matrix of normal functions given the number of spatial dimensions.

NormalInMatrixFunc    getNormalInMatrixFunc

  ( int                 rank )

{
  if ( rank == 2 )
  {
    return & get2DNormalInMatrixFuncs;
  }
  else
  {
    return & get3DNormalFuncs;
  }
}


// -----------------------------------------------------------------------
//   transformation matrix
// -----------------------------------------------------------------------  

void                  get2DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       grads,
    const Matrix&       coord )
{
  double x0 = coord(0,1) - coord(0,0);
  double y0 = coord(1,1) - coord(1,0);

  double alpha     = ::atan2 ( y0, x0 );
  double sinAlpha  = ::sin (alpha);
  double cosAlpha  = ::cos (alpha);

  Q(0,0) = - sinAlpha; Q(0,1) = cosAlpha;
  Q(1,0) =   cosAlpha; Q(1,1) = sinAlpha;
}

void                  get3DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       xgrads,
    const Matrix&       coord )
{
  using namespace jem;
  using jem::numeric::crossProduct;
  using jem::numeric::dotProduct;
  using namespace jive::geom;

  Tuple<double,3>     n, s1, s2;

  Vector              p(3);
  double              xg0, xg1;

  s1 = 0.0;
  s2 = 0.0;
  
  const int  nodeCount = coord.size (1);

  for ( int i = 0; i < nodeCount; i++ )
  {
    xg0   = xgrads(0,i);
    xg1   = xgrads(1,i);

    s1[0] += coord(0,i) * xg0;
    s1[1] += coord(1,i) * xg0;
    s1[2] += coord(2,i) * xg0;

    s2[0] += coord(0,i) * xg1;
    s2[1] += coord(1,i) * xg1;
    s2[2] += coord(2,i) * xg1;
  }

  n        = crossProduct ( s1, s2 );
  double a = ::sqrt ( dotProduct( n, n   ) );
  double b = ::sqrt ( dotProduct( s1, s1 ) );

  if ( jem::Float::isTiny( a ) || jem::Float::isTiny( b ) )
  {
    zeroVectorError ( "normal vector", "normal" );
  }

  // normalize n and s1

  n   = (1.0 / a) * n;
  s1  = (1.0 / b) * s1;

  // make s2 orthogonal to n and s1

  s2  = crossProduct ( n, s1 );

  // build the transformation matrix Q
  /*
   * Q=[
         e1.n e1.s1 e1.s2
         e2.n e2.s1 e2.s2
         e3.n e3.s1 e3.s2
       ]  

    v_loc = Q' * v_global
    I implemented Q' as Q in the following.
   */
  
  Q(0,0) = n[0]; Q(0,1) = s1[0]; Q(0,2) = s2[0];
  Q(1,0) = n[1]; Q(1,1) = s1[1]; Q(1,2) = s2[1];
  Q(2,0) = n[2]; Q(2,1) = s1[2]; Q(2,2) = s2[2]; 
  
  //System::out() << Q << "\n";
}


TransformationMatrixFunc getTransMatrixFunc

  ( int                 rank )
{
  if ( rank == 2 )
  {
    return & get2DTransMatrixFunc;
  }
  else
  {
    return & get3DTransMatrixFunc;
  }
}


//-----------------------------------------------------------------------
//   check acoustic tensor for bifurcation
//-----------------------------------------------------------------------

// Ortiz's algorithm
// tangent is the consistent material tangent which is a 4x4 matrix for plane
// strain problems.
// det(A)=det(n*CTO*n)= quartic equation in terms of tan(theta) by dividing
// det(A)=0 with cos^4(theta).

bool                      checkAcousticTensor2D 

  (       Vector&       normal,
    const Matrix&       tangent )
{
    bool localised (false);

    double D1111 = tangent(0,0);
    double D1212 = tangent(3,3);
    double D1112 = tangent(0,3);
    double D1211 = tangent(3,0);
    double D1222 = tangent(3,1);
    double D2212 = tangent(1,3);
    double D2211 = tangent(1,0);
    double D1122 = tangent(0,1);
    double D2222 = tangent(1,1);

    double a0 = D1111*D1212 - D1112*D1211; 
    double a1 = D1111*D1222 + D1111*D2212 - D1112*D2211 - D1122*D1211; 
    double a2 = D1111*D2222 + D1112*D1222 + D1211*D2212 - D1122*D1212 - D1122*D2211 - D1212*D2211; 
    double a3 = D1112*D2222 + D1211*D2222 - D1122*D2212 - D1222*D2211;
    double a4 = D1212*D2222 - D2212*D1222; 

    // finding minimum of a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
    // by finding roots of a cubic equation, which is the derivative of the
    // quartic function.

    Vector ans;
    solveCubicEqua ( ans, 4.*a4, 3.*a3, 2.*a2, a1 );

     
    const int noRoots = ans.size ( );

    double xm(0.), fm(0.);
        
    System::out() << "coefs: " << a4 << " " << a3 << " " << a2 << " "  << a1 << " " << a0 << "\n";
    System::out() << "roots: " << ans << "\n";

    for ( int i = 0; i < noRoots; i++ )
    {
        double xi = ans[i];
        double fi = a4*::pow(xi,4) + a3*::pow(xi,3) + a2*::pow(xi,2) + a1*xi + a0;
            
        System::out() << "fi: " << fi << "\n";

        // roots = [-x0,xx,x0] 
        // where -x0 and x0 both give negative fm 
        // only choose x0
        if ( ( xi > 0. ) && ( fi < fm ) )
        {
            localised = true;
            xm        = xi;
            fm        = fi;
            System::out() << "localised :)\n";
            System::out() << (180/3.14)*atan(xm) <<"\n";
        }
    }

    //xm = tan(theta);

    double theta = atan ( xm );

    normal[0] = cos ( theta );
    normal[1] = sin ( theta );

    return localised;
}

bool                      checkAcousticTensor3D 

  (       Vector&       normal,
    const Matrix&       tangent )
{
    return false;
}

CheckAcousticTensorFunc   getCheckAcousticTensorFunc

  ( int                  rank )
{
  if ( rank == 2 )
  {
    return & checkAcousticTensor2D;
  }
  else
  {
    return & checkAcousticTensor3D;
  }
}

//-----------------------------------------------------------------------
//   solveLinearEqua
//-----------------------------------------------------------------------

void                  solveLinearEqua

  ( double&       ans,
    const double  a,
    const double  b )
{
  ans = -b/a;
}

//-----------------------------------------------------------------------
//   solveQuadEqua
//-----------------------------------------------------------------------

void                  solveQuadEqua

  ( Vector&       ans,
    const double  a,
    const double  b,
    const double  c )
{
  if ( jem::numeric::abs(a) < EPS )
  {
    ans.resize( 1 );
    ans[0] = - c / b;
  }
  else
  {
    double d = b * b - 4.0 * a * c;
    
    if ( d < 0 )
    {
      using namespace jem;
      throw Error (
      JEM_FUNC,
      " imaginary principal values !!! "
      );  
    }
    else
    {
      ans.resize( 2 );
      ans[0] = 0.5 * (-b + ::sqrt( d )) / a;
      ans[1] = 0.5 * (-b - ::sqrt( d )) / a;
    }
  }
}

//-----------------------------------------------------------------------
//   solveCubicEqua
//-----------------------------------------------------------------------

void                  solveCubicEqua

  ( Vector&       ans,
    const double  a,
    const double  b,
    const double  c,
    const double  d)

{  
  // a=0: reduced to quadratic equation  

  if ( jem::numeric::abs(a) < EPS ) 
  {
    // b = 0: reduced to linear equation

    if ( jem::numeric::abs(b) < EPS )
    {
      ans.resize ( 1 );
      
      ans[0] = - d / c;

      return;
    }
    else
    {
      double D = c * c - 4.0 * b * d;
      
      if ( D < 0.0 )
      {	
	    return;
      }
      else
      {
	    D      = ::sqrt(D);

	    ans.resize ( 2 );
	    
	    ans[0] = 0.5 * (-c + D) / b;
	    ans[1] = 0.5 * (-c - D) / b;
       
	    return;
      }
    }
  }
  else
  {
    // solving cubic equation using Cardano's method.  
    // normalize to have: x^3 + aa*x^2 + bb*x + cc   
    double kk  = 1.0 / a;
    
    double aa  = b * kk;
    double bb  = c * kk;
    double cc  = d * kk;

    // change variable x=y-aa/3 to eliminate quadratic term
    // x^3 + px + q = 0

    double aa2 = aa*aa;
    double p   = ( -aa2 + 3.0 * bb ) / 9.0;
    double q   = ( 2.0 * aa * aa2 - 9.0 * aa * bb + 27.0 * cc ) / 54.0;
    double p3  = p * p * p;
    double D   = p3 + q * q;
    double aa3 = ONE_THIRD*aa;
 
    if ( jem::numeric::abs(D) < EPS )
    {
      if ( jem::numeric::abs(q) < EPS )
      {
          // one triple solution
          ans.resize(1);
          ans[0] = -aa3;
          return;
      }
      else
      {
          // one single and one double solution

          double u = pow(-q,ONE_THIRD);
          ans.resize(2);
          ans[0] = 2. * u - aa3;
          ans[1] =    - u - aa3;
          return;
      }
    }
    else
    {
        if ( D < 0. ){
            // three real solutions

            double phi = ONE_THIRD * acos ( -q / sqrt ( -p3  ) );
            double t   = 2. * sqrt ( -p );

            ans.resize(3);
            ans[0] =  t * cos ( phi       ) - aa3;
            ans[1] = -t * cos ( phi + PI3 ) - aa3;
            ans[2] = -t * cos ( phi - PI3 ) - aa3;
            return;
        }
        else{
            // one real solution

            ans.resize(1);
            double sqrtD = sqrt ( D );
            double u     = pow  ( sqrtD + jem::numeric::abs(q), ONE_THIRD );
            if ( q > 0. ){
                ans[0] = -u + p / u - aa3;
            }
            else{
                ans[0] = u - p / u - aa3;;
            }
            return;
        }
    }
  }
}


//-----------------------------------------------------------------------
//   solveCubicEqua (overloaded, Frans version)
//-----------------------------------------------------------------------

idx_t                 solveCubicEqua

  ( Tuple<double,3>& ans, 
    const double     a,
    const double     b,
    const double     c,
    const double     d)

{  
  ans = NAN;

  if ( std::abs(a) < EPS )
  {
    if ( std::abs(b) < EPS )
    {
      ans[0] = - d / c;

      return 1;
    }
    else
    {
      double D = c * c - 4.0 * b * d;
      
      if ( D < 0.0 )
      {        
        return 0;
      }
      else
      {
        D      = sqrt(D);

        ans[0] = 0.5 * (-c + D) / b;
        ans[1] = 0.5 * (-c - D) / b;
       
        return 2;
      }
    }
  }
  else
  {
    double kk  = 1.0 / a;
    
    double aa  = b * kk;
    double bb  = c * kk;
    double cc  = d * kk;

    double aa3 = aa / 3.0;

    double p, q, r;
    double phi, help;
   
    q  = ( aa * aa - 3.0 * bb ) / 9.0;
    r  = ( 2.0 * aa * aa * aa - 9.0 * aa * bb + 27.0 * cc ) / 54.0;
 
    help = r / sqrt( q * q * q );
   
    if ( std::abs (help) > 1.0 )
    {
      help = ( help < 0 ) ? -1.0 : 1.0; // prevent rounding errors
    }
   
    phi = acos ( help );
    p   = sqrt ( q    );

    ans[0] = -2.0 * p * cos ( phi / 3.0 )           - aa3;
    ans[1] = -2.0 * p * cos ( phi / 3.0 + RAD_120 ) - aa3;
    ans[2] = -2.0 * p * cos ( phi / 3.0 - RAD_120 ) - aa3;

    return 3;
  }
}


// -------------------------------------------------------
//   updateCoord2D
// -------------------------------------------------------

void            updateCoord2D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp )
{
  JEM_PRECHECK ( uCoord.size(0) == 2 &&
                 coord.size(0)  == 2 );

  uCoord(0,ALL) = coord(0,ALL) + disp[slice(0,END,2)];
  uCoord(1,ALL) = coord(1,ALL) + disp[slice(1,END,2)];
  
  //const int nnode = uCoord.size(1);
  //for (int i = 0; i < nnode; ++i )
  //{
  //  uCoord(0,i) = coord(0,i) + disp[2*i];
  //  uCoord(1,i) = coord(1,i) + disp[2*i+1];
  //}
}


// -------------------------------------------------------
//   updateCoord3D
// -------------------------------------------------------

void            updateCoord3D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp )
{
  JEM_PRECHECK ( uCoord.size(0) == 3 &&
                 coord.size(0)  == 3 );

  uCoord(0,ALL) = coord(0,ALL) + disp[slice(0,END,3)];
  uCoord(1,ALL) = coord(1,ALL) + disp[slice(1,END,3)];
  uCoord(2,ALL) = coord(2,ALL) + disp[slice(2,END,3)];
}


UpdateCoordFunc     getUpdateCoordFunc

  ( int                 rank )
{
  if     ( rank == 2 )
  {
    return & updateCoord2D;
  }
  else
  {
    return & updateCoord3D;
  }
}

//----------------------------------------------------------------------------
//   matrixFromTable
//----------------------------------------------------------------------------

void  matrixFromTable

  ( Matrix&              mat,
    const Table&         t,
    const StringVector&  colNames )

{
     Array<idx_t>    jcols    ( colNames.size() );
     Array<idx_t>    irows    ( t.rowCount() );
     idx_t           i; 

     const int rowCount = irows.size();
     const int colCount = jcols.size();
     
     mat.resize ( rowCount, colCount );
     mat = 0.0;

     for ( i = 0; i < rowCount; i++ )
     {
       irows[i] = i;
     }

     for ( i = 0; i < colCount; i++ )
     {
       jcols[i] = t.getColumnIndex ( colNames[i] );
     }

     t.findBlock ( mat, irows, jcols );

     mat.transpose();
}

void  matrixFromTable

  ( IntMatrix&     mat,
    const Table&   t,
    const String&  cols )

{
     using jem::util::StringUtils;

     StringVector colNames ( StringUtils::split( cols ) );
     Array<idx_t> jcols    ( colNames.size() );
     Array<idx_t> irows    ( t.rowCount() );
     Matrix       temp     ( irows.size(), jcols.size() );
     int          i;

     temp = 0.0;
     mat.resize ( irows.size(), jcols.size() );

     for ( i = 0; i < irows.size(); i++ )
     {
        irows[i] = i;
     }

     for ( i = 0; i < jcols.size(); i++ )
     {
        jcols[i] = t.getColumnIndex ( colNames[i] );
     }

     t.getBlock ( temp, irows, jcols );

     temp.transpose();

     //mat = jem::castTo<int> ( temp );
}

  void matrixFromTable

    ( Vector&        vec,
      const Table&   t,
      const String&  col )

{
     int  rowCount = t.rowCount();
     int  index;
     int  i;

     vec.resize ( rowCount );
     index = t.getColumnIndex ( col );

     vec = 0.0;

     for ( i = 0; i < rowCount; i++ )
     {
        vec[i] = t.getValue(i,index);
     }
}

void matrixFromTable

  ( IntVector&     vec,
    const Table&   t,
    const String&  col )

{
     int          rowCount = t.rowCount();
     int          index;
     int          i;

     vec.resize ( rowCount );
     index = t.getColumnIndex ( col );

     vec = 0;

     for ( i = 0; i < rowCount; i++ )
     {
        vec[i] = (int) t.getValue(i,index);
     }
}


// -----------------------------------------------------------------------
//   some basic functions
// -----------------------------------------------------------------------

double macaulayP ( double x ) 
{
  return 0.5*( x + std::abs (x) );
}

double macaulayN ( double x ) 
{
  return 0.5*( x - std::abs (x) );
}

double heavisideP ( double x ) 
{
  return (x < 0.) ? 0.0 : 1.0;
}

double heavisideN ( double x ) 
{
  return (x > 0.) ? 0.0 : 1.0;
}

int sign ( double x ) 
{
  return (x < 0.) ? 0 : (x > 0.);
}


// -----------------------------------------------------------------------
//   compute eigen values of 2x2 system
// -----------------------------------------------------------------------

void computeEigenValues2D 

( const Vector& stressP,
  const Vector& stressN )
{
  double sigmaxx = stressN[0];
  double sigmayy = stressN[1];
  double sigmaxy = stressN[3];

  double delta   = pow ( sigmaxx - sigmayy, 2 ) + 4. * sigmaxy * sigmaxy;

  stressP[0]     = 0.5 * ( sigmaxx + sigmayy + sqrt ( delta ) );
  stressP[1]     = 0.5 * ( sigmaxx + sigmayy - sqrt ( delta ) );
  stressP[2]     = 0.0;
}


//-----------------------------------------------------------------------
//  swapRows
//-----------------------------------------------------------------------

void  swapRows

  (  Matrix&  f,
    const int& m,
    const int& n )
  {
    int nn = f.size(1);

    double tt; //temp variable

    for ( int j = 0; j < nn; ++j ) //loop over colms
      {
        tt = f(m, j); // stote elements of m-th row 
        f(m, j) = f(n, j); // replace row-m with row-n
        f(n, j) = tt; // replace row-n with tempVar
      }

  }


//-----------------------------------------------------------------------
//  swapCols
//-----------------------------------------------------------------------

void  swapCols

  (  Matrix&  f,
    const int& m,
    const int& n )
  {
    int mm = f.size(0); 

    double tt; //temp variable

    for ( int i = 0; i < mm; ++i ) //loop over colms
      {
        tt = f(i, m); // stote elements of m-th row 
        f(i, m) = f(i, n); // replace row-m with row-n
        f(i, n) = tt; // replace row-n with tempVar
      }

  }