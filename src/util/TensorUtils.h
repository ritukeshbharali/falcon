
/** @file TensorUtils.h
 *  @brief Implements basic tensor operations.
 *  
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */

#ifndef TENSOR_UTILS_H
#define TENSOR_UTILS_H

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/Array.h>

using jem::Array;
using jem::String;
using jive::Matrix;
using jive::Vector;

typedef Array<double,4>    Quadix;

namespace tensorUtils

{

  //-----------------------------------------------------------------------
  //   public functions (common)
  //-----------------------------------------------------------------------

  double trace

    ( const Matrix&    a  );

  Matrix otimes

    ( const Vector&    a ,
      const Vector&    b ); 

  Quadix otimes

    ( const Matrix&    a ,
      const Matrix&    b );        

  Quadix otimesu

    ( const Matrix&    a ,
      const Matrix&    b );        

  Quadix otimesl

    ( const Matrix&    a ,
      const Matrix&    b );        

  double doubleDot

    ( const Matrix&    a ,
      const Matrix&    b ); 

  Matrix doubleDot

    ( const Quadix&    a ,
      const Matrix&    b ); 

  Matrix getI2

    ( const int      rnk );

  Quadix getSymI4

    ( const int      rnk );

  Quadix getSymI4Dev

    ( const int      rnk );  


  //-----------------------------------------------------------------------
  //   public functions (mechanics)
  //-----------------------------------------------------------------------  

  // 4th order tensor to 2nd order ( A_ijkl to voigt M_pq )

  Matrix tensor2voigt

    (  const Quadix&   A,
       const int       strCount );

  Vector tensor2voigtStrain

    (  const Matrix&   A,
       const int       strCount );

  Vector tensor2voigtStress

    (  const Matrix&   A,
       const int       strCount );

  Matrix voigt2tensorStrain

    (  const Vector&   c   );  

  Matrix voigt2tensorStress

    (  const Vector&   c   );

  Matrix voigt2tensorRankStrain

    (  const Vector&   c   );
    
  Matrix voigt2tensorRankStress

    (  const Vector&   c   );      

  Vector engng2tensorial

    (  const Vector&   c   );  

  Vector fill2Dstress

    (  const Vector&   a   );

  Vector fill2Dstrain

    (  const Vector&   a   ); 

  Matrix fill2Dmatrix

    (  const Matrix&   A   );      

  Vector fill3Dstress

    (  const Vector&   a   ); 

  Vector fill3Dtensorialstrain

    (  const Vector&   a   ); 

  Vector fill3Dstrain

    (  const Vector&   a   );    


}






#endif
