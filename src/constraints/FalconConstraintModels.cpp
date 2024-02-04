
#include "FalconConstraintModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------

void declareFalconConstraintModels  ()
{
  declareConstLoadArclenModel        (); 
  declareDirichletModel              ();
  declareNeumannModel                ();
  declarePeriodicModel               ();
}
