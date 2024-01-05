/** 
 *  Some basic functions common to most FE models and
 *  post-processing operations are implemented. Almost
 *  all functions are built on top of codes by written
 *  by Vinh Phu Nguyen and Frans van der Meer.
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022], 
 */


#ifndef TABLE_UTILS_H
#define TABLE_UTILS_H

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


enum ProblemType {
       PlaneStrain,
       PlaneStress,
       AxiSymmetric
};

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

#endif

