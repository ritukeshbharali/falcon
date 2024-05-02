
#include "FalconConstraintModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------

void declareFalconConstraintModels  ()
{
  declareConstLoadArclenModel        ();
  declareDispArclenModel             (); 
  declareDirichletModel              ();
  declareNeumannModel                ();
  declarePeriodicModel               ();
}
