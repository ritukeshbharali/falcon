/*
 *  Copyright (C) 2009 TU Delft. All rights reserved.
 *
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


#include <jem/util/ArrayBuffer.h>
#include <jem/base/Float.h>
#include <jem/base/array/operators.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/util/Globdat.h>
#include <jive/util/Printer.h>
#include <jive/util/Constraints.h>
#include <jive/util/XDofSpace.h>

#include "DispArclenBCs.h"
#include "util/XNames.h"

using jem::util::ArrayBuffer;
using jive::model::ActionParams;
using jive::model::Actions;

//-----------------------------------------------------------------------
//   static constants
//-----------------------------------------------------------------------

const char* DispArclenBCs::TYPE_NAME     = "DispControl";
const char* DispArclenBCs::NODES_PROP    = "nodeGroups";
const char* DispArclenBCs::MGROUP_PROP   = "imaster";
const char* DispArclenBCs::DOF_PROP      = "dofs";
const char* DispArclenBCs::LOADED_PROP   = "loaded";
const char* DispArclenBCs::LOAD_PROP     = "loads";

//-----------------------------------------------------------------------
//   constructor
//-----------------------------------------------------------------------


DispArclenBCs::DispArclenBCs ( )
{
  // initialize variables

  master_        = 0;
  loadScale_     = 0.;
  loadScale0_    = 0.;
  dispIncr_      = 0.;
  dispVal0_      = 0.;
  dispVal_       = 0.;
  loaded_        = -1;
  masterNode_    = 0;
}

DispArclenBCs::~DispArclenBCs()
{}

//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------

void DispArclenBCs::configure

  ( const Properties&  myProps,
    const Properties&  globdat )

{
  nodes_ = NodeSet::find    ( globdat );
  myProps.get( nodeGroups_, NODES_PROP );

  ngroups_ = nodeGroups_.size ( );

  // get names of dof

  myProps.get( dofTypes_, DOF_PROP );

  if ( dofTypes_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "dofTypes must have the same length as nodeGroups" );
  }

  // index of boundary with nonzero displacements

  myProps.find ( loaded_, LOADED_PROP, -1, ngroups_-1 );

  // if no such boundary specified, there must be a PointLoadModel

  if ( loaded_ < 0 ) 
  {
    master_ = -1;

    Properties loadProps;

    if ( myProps.find ( loadProps, LOAD_PROP ) )
    {
      loadBC_ = newInstance<PointLoadModel> ( LOAD_PROP );

      loadBC_-> configure( myProps, globdat );
    }
    else
    {
      throw IllegalInputException ( JEM_FUNC, 
          "no nonzero boundary condition specified" );
    }
  }
  else
  {
    loadBC_ = NIL;

    // Optionally: specify master node with additional nodegroup name.
    // This NodeGroup must have size 1 and refer to a member of the loaded
    // NodeGroup.
    // This is needed for the combination with periodic boundary condition
    // when the master node must be the corner node
    
    Assignable<NodeGroup>  group;
    String     context = "DispArclenBCs";

    group   = NodeGroup::get ( nodeGroups_[loaded_], nodes_, globdat, context );
    IdxVector iloaded ( group.getIndices () );

    if ( myProps.find ( mgroup_, MGROUP_PROP ) )
    {
      group   = NodeGroup::get ( mgroup_, nodes_, globdat, context );
      IdxVector imaster ( group.getIndices () );

      JEM_PRECHECK ( imaster.size() == 1 );
      JEM_PRECHECK ( testany ( imaster[0] == iloaded ) ); 

      masterNode_ = imaster[0];
    }
    else
    {
      masterNode_ = iloaded[0];
    }
  }
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DispArclenBCs::getConfig 

  ( const Properties& myConf,
    const Properties& globdat ) const

{
  myConf.set ( NODES_PROP,    nodeGroups_ );
  myConf.set ( DOF_PROP,      dofTypes_   );
  myConf.set ( LOADED_PROP,   loaded_     );
  myConf.set ( MGROUP_PROP,   mgroup_     );

  if ( loadBC_ != NIL )
  {
    loadBC_->getConfig ( myConf, globdat );
  }
}

//-----------------------------------------------------------------------
//   toDispControl
//-----------------------------------------------------------------------

void DispArclenBCs::toDispControl

  ( const Properties&  params,
    const Properties&  globdat )

{
  bool accepted;

  globdat.get ( accepted, "var.accepted" );

  bool temporary = false;

  params.find ( temporary, XActions::TEMPORARY );

  if ( loaded_ >= 0 )
  {
    // erase linear constraints defined during arclength control

    IdxVector  inodes;
    Assignable<NodeGroup> group;

    idx_t ig = loaded_;
    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    idx_t nn = group.size ();
    inodes . resize ( nn );
    inodes = group.getIndices ();

    idx_t itype  = dofs_->findType ( dofTypes_[ig] );

    for ( idx_t in = 0; in < nn; in++ )
    {
      if ( inodes[in] != masterNode_ )
      {
        idx_t idof = dofs_->findDofIndex ( inodes[in], itype );

        if ( idof >= 0 )
        {
          cons_->eraseConstraint ( idof );
        }
      }
    }

    if ( ! accepted )
    {
      // dispVal_ has to be set, advance() will not do anything

      if ( temporary )
      {
        // switching after crack propagation (tmpNonlin)
        // fix the current displacement 

        Vector  disp;

        StateVector::get ( disp,  dofs_, globdat );

        dispVal_ = disp[master_];
      }
      else
      {
        // switching after cancel
        // apply increment to old displacement

        dispVal_ = dispVal0_ + dispIncr_;
      }
    }
  }
  else
  {
    if ( ! ( accepted || temporary ) )
    {
      loadScale_ = loadScale0_ + dispIncr_;
    }
  }
}

//-----------------------------------------------------------------------
//   getLoadedDofs
//-----------------------------------------------------------------------

idx_t  DispArclenBCs::getLoadedDofs

  ( IdxVector&          idofs )        const

{
  // get indices of degrees of freedom with nonzero entries in 
  // external force vector (currently always a single dof)

  ArrayBuffer<idx_t> itmp;

  itmp.pushBack ( master_ );

  idofs.ref ( itmp.toArray() );

  return idofs.size();
}

//-----------------------------------------------------------------------
//   reduceDispIncr
//-----------------------------------------------------------------------

double DispArclenBCs::reduceDispIncr

  ( const double       reduction )

{
  // reduce the displacement increment and use the new increment
  // to update the prescribed value for this time step

  dispIncr_ *= reduction;

  if ( loaded_ >= 0 )
  {
    dispVal_   = dispVal0_ + dispIncr_;
  }
  else
  {
    loadScale_ = loadScale0_ + dispIncr_;
  }

  // return the new value to check for bounds in DispArclenModel

  return dispIncr_;
}


//-----------------------------------------------------------------------
//   initConstraints
//-----------------------------------------------------------------------

void DispArclenBCs::initConstraints

  ( const Properties&  globdat )

{
  idx_t                 nn;
  Assignable<NodeGroup> group;
  IdxVector             inodes;

  // Get nodes, then dofs of nodes, and constraints of dofs

  dofs_  = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat );

  String  context = nodes_.getContext();

  // loop over node groups

  for ( idx_t ig = 0; ig < ngroups_; ig++ )
  {
    group  = NodeGroup::get ( nodeGroups_[ig], nodes_, globdat, context );

    nn     = group.size();

    inodes . resize ( nn );
    inodes = group.getIndices ();

    idx_t itype  = dofs_->findType ( dofTypes_[ig] );

    if ( ig != loaded_ )  
    {
      // apply u=0 boundary conditions

      for ( idx_t in = 0; in < nn; in++ )
      {
        idx_t idof = dofs_->findDofIndex ( inodes[in], itype );

        if ( idof >= 0 )
        {
          cons_->addConstraint ( idof, 0. );
        }
      }
    }
    else
    {
      // get master dof index

      master_ = dofs_->getDofIndex ( masterNode_, itype );
    }
  }

  // compress for more efficient storage

  cons_->compress();
}

//-----------------------------------------------------------------------
//   setInitVal
//-----------------------------------------------------------------------

void DispArclenBCs::setInitVal

  ( const double initDisp )

{
  dispVal_ = dispVal0_ = initDisp;
}

//-----------------------------------------------------------------------
//   toArclControl
//-----------------------------------------------------------------------

// remove constraints already defined
// if a node group has only one node, then apply a force there
// otherwise, tie nodes together and apply a force on one 
// master node.

void DispArclenBCs::toArclControl

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::util::Printer;

  if ( loaded_ >= 0 )
  {
    IdxVector  inodes, idofs;
    Assignable<NodeGroup> group;

    IdxVector itypes(1);

    idx_t ig = loaded_;

    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    idx_t nn = group.size ();
    inodes . resize ( nn );
    idofs  . resize ( nn );
    inodes = group.getIndices ();

    idx_t itype  = itypes[0] = dofs_->findType ( dofTypes_[ig] );

    // get dofs associated with inodes
    // then erase their constraints

    idx_t ndofs = dofs_->collectDofIndices ( idofs, inodes, itypes );
    cons_->eraseConstraints ( idofs[slice(0,ndofs)] );

    master_ = dofs_->getDofIndex ( masterNode_, itype );

    // System::out() << "tying " << inodes[slice(0,END)] << " to " << masterNode_ << endl;
    // cons_->printTo(jive::util::Printer::get());
    for ( idx_t in = 0; in < nn; in++ )
    {
      if ( inodes[in] != masterNode_ )
      { 
        idx_t idof = dofs_->findDofIndex ( inodes[in], itype ) ;

        if ( idof >= 0)
        {
          cons_->addConstraint ( idof, master_, 1.0 );
        }
      }
    }
  }

  // give right initial value to ArclenSolver: 
  // ( when accepted and converged, this function is called after commit
  //   so loadScale0_ has been updated already, when switching due to 
  //   non-convergence, the old value should be used, so loadScale0_
  //   always gives the right value )

  params.set( "OldLoadScale", loadScale0_ );
}

//-----------------------------------------------------------------------
//   applyConstraints
//-----------------------------------------------------------------------

void DispArclenBCs::applyConstraints

  ( const Properties&  params,
    const Properties&  globdat )

{
  if ( loaded_ >= 0 )
  {
    // Update nonzero displacement on loaded boundary

    IdxVector             inodes;
    Assignable<NodeGroup> group;

    // the displacement increment is adjusted outside
    // precisely in the DispArclenModel

    idx_t  itype, idof;

    idx_t ig = loaded_;
    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    idx_t nn = group.size ();
    inodes . resize ( nn );
    inodes = group.getIndices ();

    itype  = dofs_->findType ( dofTypes_[ig] );

    // get amount of displacement

    for ( idx_t in = 0; in < nn; in++ )
    {
      idof = dofs_->findDofIndex ( inodes[in], itype );

      if ( idof >= 0 )
      {
        cons_->addConstraint ( idof, dispVal_ );
      }
    }

    // compress for more efficient storage

    cons_->compress();
  }
}

//-----------------------------------------------------------------------
//   storeLoadScale
//-----------------------------------------------------------------------


void DispArclenBCs::storeLoadScale

  ( const Properties&  globdat,
    const Vector&      fint )

{
  if ( loaded_ >= 0 )
  {
    // compute the total load on the loaded boundary
    // (when a PointLoadModel is defined, loadScale is updated in advance)

    Assignable<NodeGroup> group;
    IdxVector             inodes;

    idx_t ig = loaded_;
    group  = NodeGroup::find ( nodeGroups_[ig], nodes_, globdat );

    idx_t nn = group.size ();
    inodes . resize ( nn );
    inodes = group.getIndices ();

    idx_t itype = dofs_->findType ( dofTypes_[ig] );

    loadScale_ = 0.;

    for ( idx_t in = 0; in < nn; in++ )
    {
      idx_t idof = dofs_->findDofIndex ( inodes[in], itype );

      if ( idof >= 0 )
      {
        loadScale_ += fint[idof];
      }
    }
  }
}

//-----------------------------------------------------------------------
//   getUnitLoad
//-----------------------------------------------------------------------

void DispArclenBCs::getUnitLoad

  ( const Properties&  params,
    const Properties&  globdat )

{
  // get unit force vector

  if ( loadBC_ == NIL )
  {
    // trivial for single loaded node with uniform displacement

    Vector       fext;

    params.get ( fext, ActionParams::EXT_VECTOR );

    fext[master_] = 1.;
  }
  else
  {
    // otherwise, get vector defined in PointLoadModel

    params.set ( ActionParams::SCALE_FACTOR, 1. );

    loadBC_->takeAction ( Actions::GET_EXT_VECTOR, params, globdat );
  }
}

//-----------------------------------------------------------------------
//   commit
//-----------------------------------------------------------------------

// store the converged loadscale in globdat during disp control
// and the converged displacement during arclength control

void DispArclenBCs::commit

  ( const Properties&  globdat )

{
  // store converged boundary quantities 

  dispVal0_   = dispVal_;

  loadScale0_ = loadScale_;
}

//-----------------------------------------------------------------------
//   advance
//-----------------------------------------------------------------------

void DispArclenBCs::advance

// increment the displacement value for next step
// (called in displacement control at beginning of time step)

  ( const Properties&  globdat )

{
  bool accepted;

  globdat.get ( accepted, "var.accepted" );

  if ( accepted )
  {
    if ( loaded_ >= 0 )
    {
      dispVal_   = dispVal0_ + dispIncr_;
    }
    else
    {
      loadScale_ = loadScale0_ + dispIncr_;
    }
  }
}

//-----------------------------------------------------------------------
//   getLastDispIncr 
//-----------------------------------------------------------------------

double DispArclenBCs::getLastDispIncr 

  ( const Properties&  globdat ) const 

{ 
  if ( master_ >= 0 )
  {
    Vector  disp;

    StateVector::get ( disp,  dofs_, globdat );

    return ( disp[master_] - dispVal0_ );
  }
  else
  {
    return 0.;
  }
}


