/*
 *  Copyright (C) 2009 TU Delft. All rights reserved.
 *
 *  This class implements a model to define imposed displacements 
 *  on certain node groups in a computation with displacement control 
 *  and arclength control. It is used as a child of class 
 *  DispArclenModel and is able to switch between the two different 
 *  strategies. Multiple u=0 boundaries can be defined, and 1 with 
 *  prescribed nonzero displacements. For this loaded boundary, the
 *  boundary conditions are adapted when switching. Under displacement 
 *  control, displacement values are applied directly, under arclength 
 *  control, all nodes are tied to one master node on the boundary and 
 *  a force is applied to this master node. This way we can have a 
 *  unit force vector, and a load scale, which is required for the 
 *  arclength solution strategy. 
 *
 *  Authors:
 *
 *   Vinh Phu Nguyen, V.P.Nguyen@tudelft.nl
 *   Frans van der Meer, F.P.vanderMeer@tudelft.nl
 *
 *  Date:
 *   
 *   27 February 2009
 *
 */

#ifndef DISP_CONTROL_BCS_H
#define DISP_CONTROL_BCS_H


#include <jem/base/array/utilities.h>
#include <jem/base/array/select.h>
#include <jem/base/System.h>
#include <jem/util/Properties.h>
#include <jem/base/Error.h>

#include <jive/fem/NodeGroup.h>
#include <jive/model/Actions.h>
#include <jive/model/StateVector.h>
#include <jive/model/PointLoadModel.h>
#include <jive/util/Assignable.h>

namespace jive
{
  namespace util
  {
    class Constraints;
    class XDofSpace;
  }
}

using namespace jem;

using jem::Error;
using jem::util::Properties;
using jem::io::endl;
using jive::Vector;
using jive::IdxVector;
using jive::StringVector;
using jive::fem::NodeGroup;
using jive::fem::NodeSet;
using jive::model::StateVector;
using jive::model::PointLoadModel;
using jive::util::Assignable;
using jive::util::Constraints;
using jive::util::XDofSpace;

//=======================================================================
//     class DispArclenBCs
//=======================================================================

/** @brief Implements boundary constraints for arc-length models
 * 
 *  The class \c DispArclenBCs enforces a dirichlet-type arc-length model
 *  where the load on an edge is applied on one node, and all nodes on the
 *  edge are tied to that node in the context of dofs.
 *   
 *  Original Authors: Vinh Phu Nguyen    (27 February 2009, TU Delft)
 *                    Frans van der Meer (27 February 2009, TU Delft)  
 */

class DispArclenBCs : public Object
{
 public:

  typedef DispArclenBCs Self;

  static const char*         TYPE_NAME;
  static const char*         NODES_PROP;
  static const char*         MGROUP_PROP;
  static const char*         DOF_PROP;
  static const char*         LOADED_PROP;
  static const char*         LOAD_PROP;

  enum                       Control { LOAD, DISP };

                       DispArclenBCs ();

  virtual void         configure

    ( const Properties&   myProps,
      const Properties&   globdat );

  virtual void         getConfig

    ( const Properties&   myConf,
      const Properties&   globdat )      const;

  void                 toDispControl

    ( const Properties&  params,
      const Properties&  globdat );

  void                 initConstraints

    ( const Properties&   globdat );

  void                 setInitVal

    ( const double        initDisp );

  void                 toArclControl

    ( const Properties&  params,
      const Properties&  globdat );

  void                 applyConstraints

    ( const Properties&  params,
      const Properties&  globdat );

  void                 getUnitLoad

    ( const Properties&  params,
      const Properties&  globdat );

  void                 storeLoadScale

    ( const Properties&   globdat,
      const Vector&       fint );

  void                 commit

    ( const Properties&  globdat );

  void                 advance

    ( const Properties&  globdat );

  double               getLastDispIncr  
    
    ( const Properties&   globdat )      const;

  idx_t                getLoadedDofs

    ( IdxVector&          idofs )        const;

  double               reduceDispIncr

    ( const double        reduction );

  inline  void         setDispValue

    ( const double        amount );

  inline  void         setLoadScale

    ( const double        amount );

  inline  void         setDispIncr

    ( const double        val );

  inline double        getLoadScale ()  const;

  inline double        getLoadScale0()  const;

  inline double        getDispValue ()  const;

  inline double        getDispIncr  ()  const;

  void                 undoIncrement

    ( const Properties&   globdat );

 protected:

  virtual              ~DispArclenBCs ();

 private:

  Ref<PointLoadModel>     loadBC_;
  Assignable<NodeSet>     nodes_;
  Control                 control_;

  String                  mgroup_;
  Ref<XDofSpace>          dofs_;
  Ref<Constraints>        cons_;

  idx_t                   master_;
  idx_t                   masterNode_;

  double                  dispIncr_;
  double                  dispVal_;
  double                  dispVal0_;
  double                  loadScale_; 
  double                  loadScale0_;

  idx_t                   ngroups_;

  /* the following members are constant input variables 
   *
   * nodeGroups   node groups for boundary conditions
   * dofTypes     dof types of boundary condition
   * loaded       index of boundary where nonzero disp is applied
   * 
   * nodeGroups.size() == dofTypes.size()
   * 0 <= loaded < nodeGroups.size()
   */

  StringVector            nodeGroups_;
  StringVector            dofTypes_;
  idx_t                   loaded_;
};

// ======================================================
//     implementation of inline functions
// ======================================================


inline void DispArclenBCs::setDispValue

 ( const double amount )
{
  if ( loaded_ >= 0 ) dispVal_ = amount;
}

inline void DispArclenBCs::setLoadScale

 ( const double amount )
{
  if ( loadBC_ != NIL ) loadScale_ = amount;
}

inline void DispArclenBCs::setDispIncr

 ( const double val )
{
  dispIncr_ = val;
}

inline double DispArclenBCs::getLoadScale  () const 
{ 
  return loadScale_; 
}

inline double DispArclenBCs::getLoadScale0 () const 
{ 
  return loadScale0_; 
}

inline double DispArclenBCs::getDispValue  () const 
{ 
  return dispVal_;
}

inline double DispArclenBCs::getDispIncr   () const 
{ 
  return dispIncr_;  
}

inline void DispArclenBCs::undoIncrement   

  ( const Properties&  globdat ) 

{ 
  // a negative displacement increment is applied to be canceled in advance

  bool accepted;

  globdat.get ( accepted, "var.accepted" );

  if ( accepted )
  {
    if ( loaded_ >= 0 ) 
    { 
      dispVal0_ -= dispIncr_;
    }
    else
    {
      loadScale0_ -= dispIncr_;
    }
  }
}

#endif



