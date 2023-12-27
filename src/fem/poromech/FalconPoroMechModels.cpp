
#include "FalconPoroMechModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareFalconPoroMechModels                ()
{
  declareSaturatedPorousModel                   ();
  declareSaturatedPorousFractureModel           ();
  declareSaturatedPorousMicroFractureModel      ();
  declareSaturatedPorousMicroFractureExtItModel ();
  declareTwoPhaseUnsaturatedPorousModel         ();
}