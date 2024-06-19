
#ifndef FALCON_CONSTRAINT_MODELS_H
#define FALCON_CONSTRAINT_MODELS_H


//-----------------------------------------------------------------------
//   public functions
//-----------------------------------------------------------------------

/** @brief Declares all user constraint models
 *  
 *  Currently available models:
 *     Dirichlet
 *     Neumann
 *     Periodic
 */ 

void  declareFalconConstraintModels    ();

void  declareConstLoadArclenModel      ();
void  declareDispArclenModel           ();
void  declareDirichletModel            ();
void  declareNeumannModel              ();
void  declarePeriodicModel             ();

#endif


