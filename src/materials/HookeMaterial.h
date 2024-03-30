#ifndef HOOKE_MATERIAL_H
#define HOOKE_MATERIAL_H

#include "Material.h"

/*
enum ProblemType {
  PlaneStrain,
  PlaneStress,
  AxiSymmetric
};
*/

//-----------------------------------------------------------------------
//   class HookeMaterial
//-----------------------------------------------------------------------

/** @brief 
 *  The HookeMaterial class implements an isotropic linear elastic material,
 *  and is derived from the Material base class.
 */ 

class HookeMaterial : public Material
{
 public:

  static const char*      YOUNG_PROP;
  static const char*      POISSON_PROP;
  static const char*      STATE_PROP;
  static const char*      AREA_PROP;

  explicit               HookeMaterial

    ( const idx_t          rank,
      const Properties&    globdat );

  virtual void           configure

    ( const Properties&    props,
      const Properties&    globdat );

  virtual void           getConfig

    ( const Properties&    conf,
      const Properties&    globdat )   const;

  virtual void           updateConfig ();

  virtual void           update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint );

  virtual void           allocPoints

    ( const idx_t           npoints );

  virtual Matrix         getStiffMat  () const;

  virtual ProblemType    getState     () const;    

  virtual Ref<Material>  clone ( ) const;

 protected:

  virtual               ~HookeMaterial();

 protected:

  double                  area_;
  String                  stateString_;
  ProblemType             state_;

 //private:  

  Matrix                  stiffMat_;

 protected:

  void                    computeStiffMat_();

};

#endif
