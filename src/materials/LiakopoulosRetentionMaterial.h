
/** @file LiakopoulosRetentionMaterial.h
 *  @brief Implements a Liakopoulos fluid retention model.
 * 
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 20 October 2022
 * 
 *  Updates (when, what and who)
 *     - [XX YYYY 2022]
 */

#ifndef LIAKOPOULOS_RETENTION_MATERIAL_H
#define LIAKOPOULOS_RETENTION_MATERIAL_H

#include "util/BasicUtils.h"
#include "RetentionMaterial.h"

// =======================================================
//  class LiakopoulosRetentionMaterial
// =======================================================

class LiakopoulosRetentionMaterial : public RetentionMaterial
{
 public:

  static const char*      RHO_PROP;

  explicit                LiakopoulosRetentionMaterial

    ( const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const;

  virtual void           update

    ( double&              Sf,
      double&              dSfdpf,
      double&              krf,
      double&              dkrfdpf,
      const double&        pf      );

  virtual void           update

    ( double&              Sf,
      double&              dSfdpf,
      double&              krf,
      double&              dkrfdpf,
      const double&        pf,
      const double&        pg      );

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

  virtual Ref<RetentionMaterial>   clone ( ) const;    

 protected:

  virtual                ~LiakopoulosRetentionMaterial   ();

 protected:

  // Material properties

  double                  aVG_;
  double                  eVG_;
  double                  jVG_;
  double                  sRes_;

};

#endif 
