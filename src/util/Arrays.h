#ifndef ARRAYS_H
#define ARRAYS_H

#include <jem/base/Array.h>
#include <jem/base/Tuple.h>
#include <jive/Array.h>

//=======================================================================
//   using declarations
//=======================================================================

using jem::idx_t;
using jem::Array;
using jem::Tuple;
using jive::Vector;
using jive::IntVector;
using jive::Matrix;
using jive::Cubix;
using jive::IdxVector;
using jive::IntMatrix;
using jive::StringVector;

//=======================================================================
//   typedefs (Tuples)
//=======================================================================

// Some handy alias(es) to avoid some typing

typedef Tuple<double,3>       Vec3;
typedef Tuple<double,4>       Vec4;
typedef Tuple<double,6>       Vec6;
typedef Tuple<double,3,3>     Mat3;
typedef Tuple<double,4,4>     Mat4;
typedef Tuple<double,6,6>     Mat6;
typedef Tuple<double,6,3>     Mat63;

//=======================================================================
//   typedefs (Arrays)
//=======================================================================

// Some handy alias(es) to avoid some typing

typedef Array<double,4> Quadix;

#endif
