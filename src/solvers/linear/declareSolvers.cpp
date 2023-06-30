/*
 *  Copyright (C) 2019 DRG. All rights reserved.
 *
 *  This file is part of Jive, an object oriented toolkit for solving
 *  partial differential equations.
 *
 *  Commercial License Usage
 *
 *  This file may be used under the terms of a commercial license
 *  provided with the software, or under the terms contained in a written
 *  agreement between you and DRG. For more information contact DRG at
 *  http://www.dynaflow.com.
 *
 *  GNU Lesser General Public License Usage
 *
 *  Alternatively, this file may be used under the terms of the GNU
 *  Lesser General Public License version 2.1 or version 3 as published
 *  by the Free Software Foundation and appearing in the file
 *  LICENSE.LGPLv21 and LICENSE.LGPLv3 included in the packaging of this
 *  file. Please review the following information to ensure the GNU
 *  Lesser General Public License requirements will be met:
 *  https://www.gnu.org/licenses/lgpl.html and
 *  http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
 *  This file is part of Jive, an object oriented toolkit for
 *  solving partial differential equations.
 *
 *  Jive version: 3.0
 *  Date:         Fri 20 Dec 14:30:12 CET 2019
 */

/*
 * 
 *  Updates (when, what and who)
 *     - [03 May 2022]  Added external solvers Umfpack and Intel Pardiso (RB)
 *     - [30 June 2023] Added external solver library AMGCL (RB)
 */


#include <jem/base/Once.h>
#include <jive/solver/CG.h>
#include <jive/solver/GCR.h>
#include <jive/solver/GMRES.h>
#include <jive/solver/AGMRES.h>
#include <jive/solver/FGMRES.h>
#include <jive/solver/SparseLU.h>
#include <jive/solver/SkylineLU.h>
#include <jive/solver/LocalSolver.h>
#include <jive/solver/SchurSolver.h>
#include <jive/solver/SkylineSolver.h>
#include <jive/solver/VerboseSolver.h>
#include <jive/solver/DualPrecon.h>
#include <jive/solver/DummyPrecon.h>
#include <jive/solver/CoarsePrecon.h>
#include <jive/solver/SolverPrecon.h>
#include <jive/solver/NeumannPrecon.h>
#include <jive/solver/DiagPrecon.h>
#include <jive/solver/SparseILUn.h>
#include <jive/solver/SparseILUd.h>
#include <jive/solver/MultiRestrictor.h>
#include <jive/solver/SimpleRestrictor.h>
#include <jive/solver/UserdefRestrictor.h>
#include <jive/solver/NullSpaceRestrictor.h>
#include <jive/solver/RigidBodyRestrictor.h>
#include <jive/solver/declare.h>

#if defined(WITH_UMFPACK)
#include "UmfpackSolver.h"
#endif

#if defined(WITH_PARDISO)
#include "PardisoSolver.h"
#endif

#if defined(WITH_MUMPS)
#include "MUMPSSolver.h"
#endif

#if defined(WITH_AMGCL)
#include "AMGCLSolver.h"
#endif


JIVE_BEGIN_PACKAGE( solver )


//-----------------------------------------------------------------------
//   declareSolvers
//-----------------------------------------------------------------------


static void declareSolvers_ ()
{
  CG                  :: declare ();
  GCR                 :: declare ();
  GMRES               :: declare ();
  AGMRES              :: declare ();
  FGMRES              :: declare ();
  SparseLU            :: declare ();
  SkylineLU           :: declare ();
  LocalSolver         :: declare ();
  SchurSolver         :: declare ();
  SkylineSolver       :: declare ();
  VerboseSolver       :: declare ();

  DualPrecon          :: declare ();
  DummyPrecon         :: declare ();
  CoarsePrecon        :: declare ();
  SolverPrecon        :: declare ();
  NeumannPrecon       :: declare ();
  DiagPrecon          :: declare ();
  SparseILUn          :: declare ();
  SparseILUd          :: declare ();

  MultiRestrictor     :: declare ();
  SimpleRestrictor    :: declare ();
  UserdefRestrictor   :: declare ();
  NullSpaceRestrictor :: declare ();
  RigidBodyRestrictor :: declare ();

  // Added user interfaces to external solvers

  #if defined(WITH_PARDISO)
  PardisoSolver       :: declare ();
  #endif

  #if defined(WITH_UMFPACK)
  UmfpackSolver       :: declare ();
  #endif

  #if defined(WITH_MUMPS)
  MUMPSSolver         :: declare ();
  #endif

  #if defined(WITH_AMGCL)
  AMGCLSolver         :: declare ();
  #endif

}


void declareSolvers ()
{
  static jem::Once once = JEM_ONCE_INITIALIZER;

  jem::runOnce ( once, declareSolvers_ );
}


JIVE_END_PACKAGE( solver )
