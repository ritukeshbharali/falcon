/*
 * 
 *  Copyright (C) 2009 TU Delft. All rights reserved.
 *  
 *  This class implements a model that combines displacement and 
 *  arclen control. This model is designed for use in comination 
 *  with FlexArclenModule to have a flexible path following solver. 
 *  Boundary conditions are applied via a DispArclenBCs object that 
 *  is member of this class.
 *
 *  There are two different solution strategies, each with it's own 
 *  solver (NonlinModule or ArclenModule) and way of handling 
 *  boundary conditions. The FlexArclenModule decides which strategy 
 *  is applied and communicates with this model via takeAction calls. 
 *  
 *  The tasks of this model are:
 *   - let the DispArclenBCs apply correct boundary conditions
 *     (adaptive in step type and step size)
 *   - perform actions required for ArclenModule related to energy
 *     arclength method (compute dissipation etc)
 *   - pass information to FlexArclenModule for strategy decision 
 *     making
 *
 *  Authors: 
 *  V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  F.P. van der Meer, F.P.vanderMeer@tudelft.nl
 *         
 *  Date: 27 February 2009
 *
 */


#ifndef DISP_ARCLEN_MODEL_H
#define DISP_ARCLEN_MODEL_H

#include <jem/util/Flex.h>
#include <jem/io/Writer.h>
#include <jive/Array.h>
#include <jive/model/Model.h>
#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/ModelFactory.h>

#include "DispArclenBCs.h"

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

//=======================================================================
//   class DispArclenModel
//=======================================================================

/** @brief Implements a dissipation-based arc-length model
 * 
 *  The class \c DispArclenModel implements a generic dissipation-based
 *  arc-length model. It computes the dissipation-energy arc-length and
 *  unit external force in two ways, using global internal and external
 *  force vectors, and also via input from the FE model.
 * 
 *  Original Authors: Vinh Phu Nguyen    (27 February 2009, TU Delft)
 *                    Frans van der Meer (27 February 2009, TU Delft)  
 *  
 *  Contribution from RB's PhD thesis work added. (April 30, 2024)
 *  The dissipation-based arc-length as well as the unit force may be
 *  assembled by the FE model. In that case, \c DispArclenModel will not
 *  not use the global vectors to compute them. This is required for the
 *  phase-field fracture model in solid and poromechanics.
 * 
 */


class DispArclenModel : public Model
{
 public:

  typedef DispArclenModel   Self;
  typedef Model             Super;

  static const char*        TYPE_NAME;

  static const char*        CONSTRAINT_PROP;
  static const char*        OPT_ITER_PROP;
  static const char*        SWT_ITER_PROP;
  static const char*        SWT_ENER_PROP;
  static const char*        MIN_INCR_PROP;
  static const char*        MAX_INCR_PROP;
  static const char*        MAX_DISP_PROP;
  static const char*        MAX_DISS_PROP;
  static const char*        DISP_INCR_PROP;
  static const char*        INIT_DISP_PROP;
  static const char*        MIN_DISP_PROP;
  static const char*        LOAD_SCALE_PROP;
  static const char*        REDUCTION_PROP;
  static const char*        PREFER_DC_PROP;

  explicit                  DispArclenModel

    ( const String&           name  = "arclen",
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

  virtual                  ~DispArclenModel  ();


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

  void                      getStepSize_

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

  void                      connect_        ();

  void                      dofsChanged_    ();
  void                      consChanged_    ();

  double                    getReleasedEnergy_

    ( const Properties&  params,
      const Properties&  globdat );

  void                      toArclControl_

    ( const Properties&  params,
      const Properties&  globdat );

  void                      toDispControl_

    ( const Properties&  params,
      const Properties&  globdat );

 private:

  static const idx_t        U_LOAD_;
  static const char*        DELTA_STATE_;

  Ref<DispArclenBCs>        dispBC_;
  Ref<DofSpace>             dofs_;
  Ref<VectorSpace>          vspace_;

  idx_t                     updated_;
  Vector                    unitLoad_;
  IdxVector                 idofs_;

  // stored quantities

  double                    arcLength_;
  double                    lastArcl_;
  double                    energy0_;
  double                    totalDiss_;

  idx_t                     maxNIter_;

  // flags

  bool                      isDispControl_;
  bool                      triedLarge_;
  bool                      onceDown_;
  bool                      isTmpDisp_;
  bool                      hasDespaired_;
  bool                      tempStep_;

  // fancy info to the scrren

  Writer&                   out_;

  /* the following members are constant input variables 
   *
   * optIter:     optimal number of iterations (for adaptive stepping)
   * swtIter:     number of iterations for which to switch to arc-length
   * swtEnergy:   amount of energy for which to switch to arc-length
   * reduction:   factor for reduction of increments
   * minIncr:     minimum energy increment
   * maxIncr:     maximum energy increment
   * preferDisp:  switch to displacement control when hardening
   * dispIncr:    initial displacement increment
   * minDispIncr: minimum (absolute) displacement increment
   * maxDispVal:  maximum (absolute) displacement value (in disp.control)
   *
   */

  idx_t                     optIter_;
  idx_t                     swtIter_;
  double                    swtEner_;
  double                    reduction_;
  double                    minIncr_;
  double                    maxIncr_;
  bool                      preferDisp_;

  double                    initDisp_;  
  double                    dispIncr0_; 
  double                    minDispIncr_; 
  double                    maxDispVal_;
  double                    maxTotalDiss_;
};




//#######################################################################
//   Implementation
//#######################################################################


//-----------------------------------------------------------------------
//   getMaxIter
//-----------------------------------------------------------------------


inline idx_t DispArclenModel::getMaxIter () const
{
  return optIter_;
}


//-----------------------------------------------------------------------
//   get(Min|Max)Incr
//-----------------------------------------------------------------------


inline double DispArclenModel::getMinIncr () const
{
  return minIncr_;
}


inline double DispArclenModel::getMaxIncr () const
{
  return maxIncr_;
}



#endif
