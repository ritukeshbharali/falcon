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

#include "TableUtils.h"

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