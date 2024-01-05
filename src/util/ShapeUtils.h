/** 
 *  Some utility functions related to shape functions.
 *  Based on functions written by Vinh Phu Nguyen and
 *  Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */


#ifndef SHAPE_UTILS_H
#define SHAPE_UTILS_H

#include <jem/base/Array.h>
#include <jem/base/Tuple.h>
#include <jem/base/String.h>
#include <jem/base/Ref.h>

#include <jive/Array.h>
#include <jive/util/Table.h>

#include <functional>
#include <string>
#include <vector>
#include <map>

namespace jem
{
  namespace util
  {
    class Properties;
  }

  namespace io
  {
    class PrintWriter;
  }
}

using jem::Ref;
using jive::Vector;
using jive::Matrix;
using jive::Cubix;
using jive::IntVector;
using jive::IntMatrix;
using jem::Tuple;
using jem::String;
using jem::util::Properties;
using jem::io::PrintWriter;
using jive::util::Table;
using jive::StringVector;
using jem::idx_t;

//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

// A pointer to a function that computes the spatial derivatives of
// the interpolation matrix. This is the so-called B-matrix.
// it points to the corresponding function for 1D, 2D and 3D case.

typedef void        (*ShapeGradsFunc)

  ( const Matrix&       b,
    const Matrix&       g );

// a pointer to a function that computes the shape function matrix N
// it points to the corresponding function for 1D, 2D and 3D case.  

typedef void        (*ShapeFunc)

  ( const Matrix&       sfuncs,
    const Vector&       n );

typedef void        (*FShapeFunc)

  ( const Matrix&       sfuncs,
    const Vector&       n );

typedef void        (*NormalMatrixFunc)

  ( const Matrix&       normalM,
    const Vector&       normalV );

typedef void        (*NormalInMatrixFunc)

  ( const Matrix&       normalM,
    const Vector&       normalV );

typedef void        (*UpdateCoordFunc)

  ( Matrix&              uCoord,
    const Matrix&        coord, 
    const Vector&        disp );

typedef void       (*TransformationMatrixFunc)

  ( Matrix&             Q,
    const Matrix&       grads,
    const Matrix&       coords );


typedef void       (*DoubleShearComponentsFunc)

  ( Vector&            v );

typedef bool       (*CheckAcousticTensorFunc)

  (       Vector&        n,
    const Matrix&        tangent  );  

typedef void       (*GetIandPFunc)
  
  ( const Matrix&       P,
    const Vector&       I );


//***********************************************************************
//   public functions
//***********************************************************************

//-----------------------------------------------------------------------
//   B matrix
//-----------------------------------------------------------------------

void                  get1DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

void                  get2DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

void                  get3DShapeGrads

  ( const Matrix&       b,
    const Vector&       g );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeGradsFunc        getShapeGradsFunc

  ( int                 rank );


// -----------------------------------------------------------------------
//   N matrix
// -----------------------------------------------------------------------

void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get2DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

void                  get3DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );


// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

ShapeFunc             getShapeFunc

  ( int                 rank );


// -----------------------------------------------------------------------
//   matrix of normal vectors
// -----------------------------------------------------------------------
  
 void                  get2DNormalFuncs

  ( const Matrix&       normalM,
    const Vector&       normalV );
  
 void                  get3DNormalFuncs

  ( const Matrix&       normalM,
    const Vector&       normalV );

// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

NormalMatrixFunc        getNormalFunc

  ( int                 rank );

// -----------------------------------------------------------------------
//   matrix of normal vectors (new version 2x4 matrix for 2D)
// -----------------------------------------------------------------------
  
 void                  get2DNormalInMatrixFuncs

  ( const Matrix&       normalM,
    const Vector&       normalV );
  
// A function that returns a pointer to a function that computes the
// B-matrix given the number of spatial dimensions.

NormalInMatrixFunc      getNormalInMatrixFunc

  ( int                 rank );


// -----------------------------------------------------------------------
//   transformation matrix
// -----------------------------------------------------------------------  

void                  get2DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       grads,
    const Matrix&       coords );

void                  get3DTransMatrixFunc

  (       Matrix&       Q,
    const Matrix&       grads,
    const Matrix&       coords );

TransformationMatrixFunc getTransMatrixFunc

  ( int                 rank );


// -----------------------------------------------------------------------
//   check acoustic tensor for bifurcation
// -----------------------------------------------------------------------

bool                  checkAcousticTensor2D 

  (       Vector&       normal,
    const Matrix&       tangent );

bool                  checkAcousticTensor3D 

  (       Vector&       normal,
    const Matrix&       tangent );

CheckAcousticTensorFunc          getCheckAcousticTensorFunc

  ( int                  rank );
  

// -----------------------------------------------------------------------
//   solve linear, quadratic and cubic equations
// -----------------------------------------------------------------------

// Solving linear equation of one variable
// a*x + b = 0

void                  solveLinearEqua

  ( double&       ans,
    const double a,
    const double b );

// Solving quadratic equation
// a*x^2 + b*x + c = 0

void                  solveQuadEqua

  ( Vector&       ans,
    const double a,
    const double b,
    const double c );

// Solving cubic equation
// a*x^3 + b*x^2 + c*x + c = 0

void                  solveCubicEqua

  ( Vector&       ans,
    const double a,
    const double b,
    const double c,
    const double d);

// Solving cubic equation (overloaded, Frans version)
// a*x^3 + b*x^2 + c*x + d = 0
// only for functions with real roots
// returns the number of roots

idx_t                 solveCubicEqua

  ( Tuple<double,3>& ans,
    const double     a,
    const double     b,
    const double     c,
    const double     d); 


// -----------------------------------------------------------------------
//   update coord with displacements (deformed mesh)
// -----------------------------------------------------------------------

void            updateCoord2D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp );

void            updateCoord3D 

  ( Matrix&       uCoord,
    const Matrix& coord, 
    const Vector& disp );


UpdateCoordFunc     getUpdateCoordFunc

  ( int                 rank );


// -----------------------------------------------------------------------
//   create matrix from table
// -----------------------------------------------------------------------  

void  matrixFromTable

  ( Matrix&                       mat,
    const Table&                  t,
    const StringVector&           colNames );

void  matrixFromTable

  ( IntMatrix&                    mat,
    const Table&                  t,
    const String&                 cols );

void matrixFromTable

  ( Vector&                       vec,
    const Table&                  t,
    const String&                 col );

void matrixFromTable

  ( IntVector&                    vec,
    const Table&                  t,
    const String&                 col );


// -----------------------------------------------------------------------
//   some basic functions
// -----------------------------------------------------------------------

double macaulayP  ( double x );
double macaulayN  ( double x );
double heavisideP ( double x );
double heavisideN ( double x );
int    sign       ( double x );


// -----------------------------------------------------------------------
//   compute eigen values of 2x2 system
// -----------------------------------------------------------------------

void   computeEigenValues2D 

  ( const Vector& stressP,
    const Vector& stressN );


//-----------------------------------------------------------------------
//  swapRows
//-----------------------------------------------------------------------
void  swapRows

  (  Matrix&  f,
    const int& m,
    const int& n );


//-----------------------------------------------------------------------
//  swapCols
//-----------------------------------------------------------------------
void  swapCols

  (  Matrix&  f,
    const int& m,
    const int& n );





#endif

