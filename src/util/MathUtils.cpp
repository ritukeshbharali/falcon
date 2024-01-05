/** 
 *  Some basic functions common to most FE models and
 *  post-processing operations are implemented. Almost
 *  all functions are built on top of codes by written
 *  by Vinh Phu Nguyen and Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */

#include <cmath>

#include <jem/base/Error.h>
#include <jem/numeric/utilities.h>

#include "Constants.h"
#include "MathUtils.h"

//========================================================================
//   Basic functions
//========================================================================

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


//========================================================================
//   Solve algebraic equations
//========================================================================

//-----------------------------------------------------------------------
//   solveLinearEqua
//-----------------------------------------------------------------------

void                  solveLinearEqn

  ( double&       ans,
    const double  a,
    const double  b )
{
  ans = -b/a;
}

//-----------------------------------------------------------------------
//   solveQuadEqua
//-----------------------------------------------------------------------

void                  solveQuadEqn

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
      ans[0] = 0.5 * (-b + sqrt( d )) / a;
      ans[1] = 0.5 * (-b - sqrt( d )) / a;
    }
  }
}

//-----------------------------------------------------------------------
//   solveCubicEqua
//-----------------------------------------------------------------------

// Solving cubic equation
// a*x^3 + b*x^2 + c*x + d = 0
// only for functions with real roots
// returns the number of roots

idx_t                 solveCubicEqn

  ( Tuple<double,3>& ans, 
    const double     a,
    const double     b,
    const double     c,
    const double     d)

{  
  ans = NAN;

  if ( jem::numeric::abs(a) < EPS )
  {
    if ( jem::numeric::abs(b) < EPS )
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
   
    if ( jem::numeric::abs (help) > 1.0 )
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