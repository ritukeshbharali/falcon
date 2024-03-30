
#include "FalconSolidMechModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareFalconSolidMechModels      ()
{
  declareLinearElasticityModel         ();
  declarePhaseFractureModel            ();
  declarePhaseFractureExtModel         ();
  declarePhaseFractureExtItModel       ();
  declareMicroPhaseFractureModel       ();
  declareMicroPhaseFractureExtModel    ();
  declareMicroPhaseFractureExtItModel  ();
  declareGradientDamageModel           ();
}
