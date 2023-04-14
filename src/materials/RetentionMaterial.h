#ifndef RETENTION_MATERIAL_H
#define RETENTION_MATERIAL_H

#include <jem/base/Object.h>
#include <jive/Array.h>

using jem::Ref;
using jem::String;
using jive::Vector;
using jive::Matrix;
using jem::Object;
using jem::util::Properties;


//-----------------------------------------------------------------------
//   class RetentionMaterial
//-----------------------------------------------------------------------

/** @brief 
 * 
 *  The RetentionMaterial class is the base class, providing common 
 *  functionalities for different water retention curves. 
 * 
 *  Derived classes are expected to provide implementation for the 'update'
 *  and 'clone' functions. RetentionMaterial class is derived from base 
 *  class Object, which offers garbage collection. 
 */ 

class RetentionMaterial : public Object
{
 public:

  explicit               RetentionMaterial

    ( const Properties&    globdat );

  virtual void           configure

    ( const Properties&    props,
      const Properties&    globdat );

  virtual void           getConfig

    ( const Properties&    conf,
      const Properties&    globdat ) const;

  virtual void           update

    ( double&              Sf,
      double&              dSfdpf,
      double&              krf,
      double&              dkrfdpf,
      const double&        pf      ) = 0;

  virtual void           update

    ( double&              Sf,
      double&              dSfdpf,
      double&              krf,
      double&              dkrfdpf,
      const double&        pf,
      const double&        pg      ) = 0;  

  virtual void           updateSaturation

    ( double&              Sf,
      double&              dSfdpf,
      const double&        pf      );

  virtual void           updateSaturation

    ( double&              Sf,
      double&              dSfdpf,
      const double&        pf,
      const double&        pg      );      

  virtual void           updatePermeability

    ( double&              krf,
      double&              dkrfdpf,
      const double&        pf      );

  virtual void           updatePermeability

    ( double&              krf,
      double&              dkrfdpf,
      const double&        pf,
      const double&        pg      );        

  virtual Ref<RetentionMaterial>   clone ( ) const = 0;

  //  Getters (provides read-only/copy access private/protected members of the class)

  inline double           getRho     () const;

  //  Setters (modify private/protected members of the class)  

  inline void             setRho   ( double rho   );

 protected:

                       ~RetentionMaterial ();

    double              rho_;
    double              g_;                       
};

//=======================================================================
//   GETTERS
//=======================================================================

//-----------------------------------------------------------------------
//   getYoung
//-----------------------------------------------------------------------

inline double RetentionMaterial::getRho () const
{
  return rho_;
}

//=======================================================================
//   SETTERS
//=======================================================================

inline void RetentionMaterial::setRho ( double rho ) { rho_ = rho; }

//-----------------------------------------------------------------------
//   newRetentionMaterial
//-----------------------------------------------------------------------

Ref<RetentionMaterial>  newRetentionMaterial

    ( const String&       name,
      const Properties&   conf,
      const Properties&   props,
      const Properties&   globdat );

#endif
