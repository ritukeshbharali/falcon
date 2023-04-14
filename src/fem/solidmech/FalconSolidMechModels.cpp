
#include "FalconSolidMechModels.h"


//-----------------------------------------------------------------------
//   declareModels
//-----------------------------------------------------------------------


void declareFalconSolidMechModels      ()
{
  
  declareLinearElasticityModel         ();
  declarePhaseFractureModel            ();
  declarePhaseFractureExtModel         ();
  declareMicroPhaseFractureModel       ();
  declareMicroPhaseFractureExtModel    ();

}
