#ifndef MATERIAL_H
#define MATERIAL_H

#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/base/Object.h>
#include <jive/Array.h>
#include <jive/util/XTable.h>

using jem::System;
using jem::idx_t;
using jem::Ref;
using jem::Object;
using jem::String;
using jem::util::Properties;
using jive::Vector;
using jive::Matrix;
using jive::IdxVector;
using jive::BoolVector;
using jive::StringVector;


//-----------------------------------------------------------------------
//   class Material
//-----------------------------------------------------------------------

/** @brief 
 * 
 *  The Material class is the base class, providing common functionalities
 *  for modelling different types of materials. 
 * 
 *  Derived classes are expected to provide implementation for the 'update',
 *  'updateConfig' and 'clone' functions. Material class is derived from 
 *  base class Object, which offers garbage collection. 
 */ 

class Material : public Object
{
 public:
    explicit              Material

    ( const idx_t          rank,
      const Properties&    globdat );

  virtual void           configure

    ( const Properties&    props,
      const Properties&    globdat );

  virtual void           getConfig

    ( const Properties&    conf,
      const Properties&    globdat ) const;

  virtual void           updateConfig () = 0;  

  virtual void           update

    ( Vector&              stress,
      Matrix&              stiff,
      const Vector&        strain,
      const idx_t          ipoint ) = 0;

  virtual void           allocPoints

    ( const idx_t          npoints );

  virtual void           deallocPoints

    ( const idx_t          npoints );  

  virtual void           commit ();

  virtual void           checkCommit

    ( const Properties&     params  );

  virtual void           commitOne

    ( const idx_t           ipoint );

  virtual void           cancel ();

  virtual idx_t          pointCount () const;

  virtual bool           isLoading

    ( const idx_t           ipoint ) const;

  virtual bool           wasLoading

    ( const idx_t           ipoint ) const;

  bool                   despair ();

  void                   endDespair ();

  virtual Ref<Material>  clone ( ) const = 0;

  //  Getters (provides read-only/copy access private/protected members of the class)

  virtual void           getHistory

    ( Vector&              hvals,
      const idx_t          mpoint );

  inline void            getHistoryNames   
    
    ( const StringVector&   hnames ) const;  

  inline idx_t            getHistoryCount () const;
  inline StringVector     getHistoryNames () const;
  inline double           getYoung        () const;
  inline double           getPoisson      () const;

  //  Setters (modify private/protected members of the class)  

  virtual void           setHistory

    ( const Vector&        hvals,
      const idx_t          mpoint );

  inline void            setYoung   ( double young   );
  inline void            setPoisson ( double poisson );  

 protected:

                       ~Material ();
      
 protected:                         

  int                  rank_;
  double               young_;
  double               poisson_;

  bool                 desperateMode_;
  BoolVector           hasSwitched_;
  BoolVector           useSecant_;

  StringVector         historyNames_;
  
};

//=======================================================================
//   GETTERS
//=======================================================================

//-----------------------------------------------------------------------
//   getHistoryCount
//-----------------------------------------------------------------------

inline idx_t Material::getHistoryCount () const
{
  return historyNames_.size();
}

//-----------------------------------------------------------------------
//   getHistoryNames
//-----------------------------------------------------------------------

void Material::getHistoryNames 

  ( const StringVector&  hnames ) const
{
  JEM_ASSERT ( hnames.size() == getHistoryCount() );
  hnames = historyNames_;
}

//-----------------------------------------------------------------------
//   getHistoryNames
//-----------------------------------------------------------------------

inline StringVector Material::getHistoryNames 

  () const
{
  return historyNames_;
}

//-----------------------------------------------------------------------
//   getYoung
//-----------------------------------------------------------------------

inline double Material::getYoung  () const
{
  return young_;
}

//-----------------------------------------------------------------------
//   getPoisson
//-----------------------------------------------------------------------

inline double Material::getPoisson() const
{
  return poisson_;
}

//=======================================================================
//   SETTERS
//=======================================================================

inline void Material::setYoung   ( double young   ) { young_   = young;   }
inline void Material::setPoisson ( double poisson ) { poisson_ = poisson; }

//-----------------------------------------------------------------------
//   newMaterial
//-----------------------------------------------------------------------

Ref<Material>  newMaterial

    ( const String&       name,
      const Properties&   conf,
      const Properties&   props,
      const Properties&   globdat );

#endif
