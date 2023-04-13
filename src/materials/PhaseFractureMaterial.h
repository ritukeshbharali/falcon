/*  
 *  This pure virtual class implements functions for a phase-field
 *  fracture material model. Any phase-field fracture material 
 *  should be derived from Material as well as from Softening.
 */

#ifndef PHASE_FRACTURE_MATERIAL_H 
#define PHASE_FRACTURE_MATERIAL_H

class PhaseFractureMaterial 
{
 public:
  
  virtual double          givePsi           () const = 0;
  virtual Vector          giveDerivPsi      () const = 0;

  /**
   * @brief      { Stress update for GP 'ip' given total strain }
   *
   * @param      stress  The stress vector
   * @param      stiff   The tangent stiffness
   * @param[in]  strain  The total strain vector
   * @param[in]  gphi    The material degradation 
   * @param[in]  ip      The ID of Gauss point 
   */

  virtual void            update

    ( Vector&         stress,
      Matrix&         stiff,
      const Vector&   strain,
      const double&   gphi,
      int             ip )       = 0;

  /**
   * @brief      { Stress update for GP 'ip' given total strain. 
   *               Required for micromorphic phase-field models.}
   *
   * @param      stressP  The positive stress vector
   * @param      stressN  The negative stress vector
   * @param      stiffP   The positive tangent stiffness
   * @param      stiffN   The negative tangent stiffness
   * @param[in]  strain   The total strain vector
   * @param[in]  ip       The ID of Gauss point 
   */  

  virtual void            update

    ( Vector&         stressP,
      Vector&         stressN,
      Matrix&         stiffP,
      Matrix&         stiffN,
      const Vector&   strain,
      int             ip )       = 0;

};

#endif
