/** 
 *  @file ShapeUtils.h
 *  @brief Utility functions related to shape functions, strain-displacement 
 *  matrices, and finite strain kinematics.
 *  Based on functions written by Vinh Phu Nguyen and Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [05 October 2026] Removed functions not needed. Added functions
 *       for large deformation/finite strain. (RB)
 */


#ifndef SHAPE_UTILS_H
#define SHAPE_UTILS_H

#include <jem/base/Array.h>
#include <jem/base/Tuple.h>
#include <jem/base/String.h>
#include <jem/base/Ref.h>
#include <jem/numeric/algebra/matmul.h>
#include <jem/numeric/algebra/MatmulChain.h>

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

using jem::numeric::matmul;
using jem::numeric::MatmulChain;

typedef MatmulChain<double,1>  MChain1;
typedef MatmulChain<double,2>  MChain2;
typedef MatmulChain<double,3>  MChain3;
typedef MatmulChain<double,4>  MChain4;

//-----------------------------------------------------------------------
//   typedefs
//-----------------------------------------------------------------------

/**
 * @brief Pointer to a function computing spatial derivatives of the 
 * interpolation matrix (B-matrix).
 */
typedef void        (*ShapeGradsFunc)

  ( const Matrix&       b,
    const Matrix&       g );

/**
 * @brief Pointer to a function computing the linear finite strain
 * B-matrix.
 */
typedef void        (*BMatrixLinFunc)

  (       Matrix&       b,
    const Matrix&       f,
    const Matrix&       g );

/**
 * @brief Pointer to a function computing the interpolation or shape
 * function matrix (N-matrix).
 */
typedef void        (*ShapeFunc)

  ( const Matrix&       sfuncs,
    const Vector&       n );

//=======================================================================
//   B matrix
//=======================================================================

/**
 * @brief Computes the 1D shape function gradients (B-matrix).
 * @param b Output strain-displacement matrix (1 x n).
 * @param g Input shape function gradients (1 x n).
 */
void                  get1DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

/**
 * @brief Computes the 2D strain-displacement matrix (B-matrix).
 * @param b Output strain-displacement matrix (4 x 2n).
 * @param g Input shape function gradients (2 x n).
 */
void                  get2DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

/**
 * @brief Computes the 3D strain-displacement matrix (B-matrix).
 * @param b Output strain-displacement matrix (6 x 3n).
 * @param g Input shape function gradients.
 */
void                  get3DShapeGrads

  ( const Matrix&       b,
    const Matrix&       g );

/**
 * @brief Returns a function pointer to the appropriate shape gradients 
 * function based on spatial dimension.
 * @param rank Spatial dimension (1, 2, or 3).
 * @return ShapeGradsFunc Pointer to the corresponding 1D, 2D, or 3D 
 * shape gradients function.
 */
ShapeGradsFunc        getShapeGradsFunc

  ( int                 rank );


//=======================================================================
//   B matrix (Linear): used for Finite Strain problems
//=======================================================================

/**
 * @brief Computes the linear finite strain B-matrix in 2D.
 * @param b Output B-matrix (4 x 2n).
 * @param f Deformation gradient matrix (2 x 2).
 * @param g Shape function gradients (2 x n).
 */
void  getBMatrixLin2D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g );

/**
 * @brief Computes the linear finite strain B-matrix in 3D.
 * @param b Output B-matrix (6 x 3n).
 * @param f Deformation gradient matrix (3 x 3).
 * @param g Shape function gradients (3 x n).
 */
void  getBMatrixLin3D

  (       Matrix&  b,
    const Matrix&  f,
    const Matrix&  g );

/**
 * @brief Returns a function pointer to the appropriate linear finite 
 * strain B-matrix function.
 * @param rank Spatial dimension (2 or 3).
 * @return BMatrixLinFunc Pointer to the linear B-matrix function.
 */
BMatrixLinFunc        getBMatrixLinFunc

  ( int                 rank );

//=======================================================================
//   N matrix
//=======================================================================

/**
 * @brief Computes the 1D shape function matrix (N-matrix).
 * @param sfuncs Output shape function matrix.
 * @param n Input shape function values vector.
 */
void                  get1DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

/**
 * @brief Computes the 2D shape function matrix (N-matrix).
 * @param sfuncs Output shape function matrix (2 x 2n).
 * @param n Input shape function values vector (n).
 */
void                  get2DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );

/**
 * @brief Computes the 3D shape function matrix (N-matrix).
 * @param sfuncs Output shape function matrix (3 x 3n).
 * @param n Input shape function values vector (n).
 */
void                  get3DShapeFuncs

  ( const Matrix&       sfuncs,
    const Vector&       n );


/**
 * @brief Returns a function pointer to the appropriate shape function matrix (N-matrix) builder.
 * @param rank Spatial dimension (1, 2, or 3).
 * @return ShapeFunc Pointer to the shape function matrix builder.
 */
ShapeFunc             getShapeFunc

  ( int                 rank );


//=======================================================================
//   Finite Strain utility functions
//=======================================================================

/**
 * @brief Adds the geometric nonlinear contribution to the element stiffness 
 * matrix using Voigt stress vector.
 * @param k Element stiffness matrix.
 * @param pk2 Second Piola-Kirchhoff stress in Voigt notation.
 * @param g Shape function gradients.
 * @param w Integration weight factor.
 */
void  addGeomNonlinElemMat

  (       Matrix&  k,
    const Vector&  pk2,
    const Matrix&  g,
    const double   w );

/**
 * @brief Adds the geometric nonlinear contribution to the element stiffness 
 * matrix using tensor stress.
 * @param k Element stiffness matrix.
 * @param pk2 Second Piola-Kirchhoff stress tensor.
 * @param g Shape function gradients.
 * @param w Integration weight factor.
 */
void  addGeomNonlinElemMat

  (       Matrix&  k,
    const Matrix&  pk2,
    const Matrix&  g,
    const double   w );

/**
 * @brief Computes the deformation gradient tensor from displacement vector
 * and shape gradients.
 * @param f Output deformation gradient tensor (rank x rank).
 * @param u Nodal displacement vector.
 * @param g Shape function gradients.
 */
void                  getDeformationGradient

  (       Matrix&       f,
    const Vector&       u,
    const Matrix&       g );

/**
 * @brief Computes the Green-Lagrange strain tensor from the deformation 
 * gradient.
 * @param e Output Green-Lagrange strain tensor (rank x rank).
 * @param f Deformation gradient tensor.
 */
void                  getGreenLagrangeStrain

  (       Matrix&       e,
    const Matrix&       f );

/**
 * @brief Computes the Green-Lagrange strain in Voigt notation from the 
 * deformation gradient.
 * @param eps Output Green-Lagrange strain vector in Voigt notation.
 * @param f Deformation gradient tensor.
 */
void                  getGreenLagrangeStrain

  (       Vector&       eps,
    const Matrix&       f );

/**
 * @brief Computes the First Piola-Kirchhoff stress tensor.
 * @param pk1 Output First Piola-Kirchhoff stress tensor.
 * @param f Deformation gradient tensor.
 * @param pk2 Second Piola-Kirchhoff stress tensor.
 */
void                  getFirstPiolaStress

  (       Matrix&       pk1,
    const Matrix&       f,
    const Matrix&       pk2 );

/**
 * @brief Computes the First Piola-Kirchhoff stress tensor.
 * @param pk1 Output First Piola-Kirchhoff stress tensor.
 * @param f Deformation gradient tensor.
 * @param pk2 Second Piola-Kirchhoff stress in Voigt notation.
 */
void                  getFirstPiolaStress

  (       Matrix&       pk1,
    const Matrix&       f,
    const Vector&       pk2 );

/**
 * @brief Computes the Cauchy stress vector in Voigt notation.
 * @param sigma Output Cauchy stress vector in Voigt notation.
 * @param f Deformation gradient tensor.
 * @param pk2 Second Piola-Kirchhoff stress tensor.
 */
void                  getCauchyStress

  (       Vector&       sigma,
    const Matrix&       f,
    const Matrix&       pk2 );

/**
 * @brief Computes the Cauchy stress vector in Voigt notation.
 * @param sigma Output Cauchy stress vector in Voigt notation.
 * @param f Deformation gradient tensor.
 * @param pk2 Second Piola-Kirchhoff stress in Voigt notation.
 */
void                  getCauchyStress

  (       Vector&       sigma,
    const Matrix&       f,
    const Vector&       pk2 );

#endif

