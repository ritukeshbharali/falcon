/** @file ConstLoadArclenModel.h
 *  @brief Implements constant load arc-length model in
 *  conjunction with the TSArclenModule.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 11 January 2024
 *  
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2024] 
 */


#ifndef CONSTLOAD_ARCLEN_MODEL_H
#define CONSTLOAD_ARCLEN_MODEL_H

#include <jem/util/Flex.h>
#include <jem/io/Writer.h>
#include <jive/Array.h>
#include <jive/model/Model.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>

namespace jive
{
  namespace util
  {
    class Constraints;
    class XDofSpace;
  }
}

namespace jive
{
  namespace model
  {
    class Model;
  }
}

namespace jive
{
  namespace algebra
  {
    class VectorSpace;
  }
}

using namespace jem;

using jem::util::Properties;
using jem::util::Flex;
using jem::io::Writer;
using jive::Vector;
using jive::IdxVector;
using jive::algebra::VectorSpace;
using jive::model::Model;
using jive::StringVector;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::model::Model;
using jive::util::XDofSpace;
using jive::util::DofSpace;
using jive::util::Constraints;

//-----------------------------------------------------------------------
//   class ConstLoadArclenModel
//-----------------------------------------------------------------------

/** @brief 
 *  The ConstLoadArclenModel class operates with the time-step computing
 *  TSArclenModule. It is used for problems where scaling of the external
 *  force is not possible.
 * 
 *  Usage:
 * 
 *  arcmodel =
 *  {
 *    type      = "ConstLoadArclen";
 *    optIter   = 10;
 *    maxIncr   = 0.1;
 *    minIncr   = 1.e-4;
 *    swtEnergy = 1.e-7;
 *    swtIter   = 5;
 *  };
 */ 


class ConstLoadArclenModel : public Model
{
 public:

  typedef ConstLoadArclenModel   Self;
  typedef Model                  Super;

  static const char*        TYPE_NAME;

  static const char*        OPT_ITER_PROP;
  static const char*        SWT_ITER_PROP;
  static const char*        SWT_ENER_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        REDUCTION_PROP;

  explicit                  ConstLoadArclenModel

    ( const String&           name  = "constArclen",
      const Ref<Model>&       child = NIL );

  virtual Model*            findModel

    ( const String&           name )         const;

  virtual void              configure

    ( const Properties&       props,
      const Properties&       globdat );

  virtual void              getConfig

    ( const Properties&       conf,
      const Properties&       globdat )      const;

  virtual bool              takeAction

    ( const String&           action,
      const Properties&       params,
      const Properties&       globdat );

  void                      setMaxIter

    ( idx_t                   count );

  inline idx_t              getMaxIter      () const;

  void                      setLoadIncr

    ( double                  incr );

  inline double             getLoadIncr     () const;

  void                      setIncrRange

    ( double                  minIncr,
      double                  maxIncr );

  inline double             getMinIncr      () const;
  inline double             getMaxIncr      () const;


  static Ref<Model>         makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       globdat );


 protected:

  virtual                  ~ConstLoadArclenModel  ();


 private:

  void                      init_

    ( const Properties&       globdat );

  void                      initLoad_

    ( const Properties&       globdat );

  void                      evalArcFunc_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getUnitLoad_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      getExtVector_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      commit_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      cancel_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      reportProgress_  ()  const;

  void                      checkCommit_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      checkSwitch_

    ( const Properties&       params,
      const Properties&       globdat );

  double                    getStepSize_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      setStepSize_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      reduceStep_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      increaseStep_

    ( const Properties&       params,
      const Properties&       globdat );

  void                      dofsChanged_    ();
  void                      consChanged_    ();

  void                      toArclControl_

    ( const Properties&  params,
      const Properties&  globdat );

  void                      toLoadControl_

    ( const Properties&  params,
      const Properties&  globdat );

 private:

  static const idx_t        U_LOAD_;
  static const char*        DELTA_STATE_;

  Ref<DofSpace>             dofs_;
  idx_t                     updated_;

  // input params

  idx_t                     optIter_;
  idx_t                     swtIter_;
  double                    swtEner_;
  double                    reduction_;
  double                    minIncr_;
  double                    maxIncr_;

  // arc-length constraint related quantities

  double                    arcLength_;
  double                    lastArclength_;

  double                    arcLength0_;
  double                    gamma_;

  // time step

  double                    dtime_;
  double                    dtime0_;

  // flags

  bool                      isLoadControl_;
  bool                      onceDown_;
  bool                      isTmpLoad_;
  bool                      hasDespaired_;
  bool                      tempStep_;
  bool                      triedLarge_;

  // statistics

  idx_t                     maxNIter_;

  // fancy info to the screen

  Writer&                   out_;

};



//#######################################################################
//   Implementation
//#######################################################################


//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


inline idx_t ConstLoadArclenModel::getMaxIter () const
{
  return optIter_;
}


//-----------------------------------------------------------------------
//   get(Min|Max)Arc-length
//-----------------------------------------------------------------------


inline double ConstLoadArclenModel::getMinIncr () const
{
  return minIncr_;
}


inline double ConstLoadArclenModel::getMaxIncr () const
{
  return maxIncr_;
}


#endif
