
#include "FalconPoroMechModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareFalconPoroMechModels               ()
{
	declareSaturatedPorousModel                  ();
	declareSaturatedPorousFractureModel          ();
  declareSaturatedPorousMicroFractureModel     ();
	declareTwoPhaseUnsaturatedPorousModel        ();
}