
/** @file PeriodicModel.cpp
 *  @brief Implements Periodic boundary conditions.
 *  
 *  This class implements a model for setting Periodic 
 *  boundary conditions for certain group of nodes. It
 *  must be included after the FE Model in *.pro file.
 *  The step size set in this model may be modified by
 *  any module that throws the action 
 *  SolverNames::SET_STEP_SIZE.
 * 
 *  Equation: nodeGroup1 = nodeGroup0 + scale * stepsize
 *  Jive:       idof        jdof      + coeff * rval            
 */

/* Include jem and jive headers */

#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/model/Actions.h>

/* Include user-defined class headers */

#include "PeriodicModel.h"
#include "util/XNames.h"

using jive::model::Actions;


//=======================================================================
//   class PeriodicModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  PeriodicModel::TYPE_NAME          = "Periodic";
const char*  PeriodicModel::NODE_GROUPS0_PROP  = "nodeGroups0";
const char*  PeriodicModel::NODE_GROUPS1_PROP  = "nodeGroups1";
const char*  PeriodicModel::DOFS_PROP          = "dofs";
const char*  PeriodicModel::FACTORS_PROP       = "factors";
const char*  PeriodicModel::STEP_SIZE_PROP     = "stepSize";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


PeriodicModel::PeriodicModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  Properties    myProps = props.findProps ( myName_ );
  Properties    myConf  = props.makeProps ( myName_ );

  const String  context = getContext ();

  // Get node groups assigned to this model.

  myProps.get( nodeGroups0_, NODE_GROUPS0_PROP );
  myConf. set( NODE_GROUPS0_PROP, nodeGroups0_ );

  myProps.get( nodeGroups1_, NODE_GROUPS1_PROP );
  myConf. set( NODE_GROUPS1_PROP, nodeGroups1_ );

  if ( nodeGroups0_.size() != nodeGroups1_.size() )
  {
    throw IllegalInputException ( JEM_FUNC,
          "nodeGroups0 must have the same length as nodeGroups1" );
  }

  ngroups_ = nodeGroups0_.size();

  /*
   *  Get dofs associated with the node groups. The length of
   *  dofs must be equal to the length of node groups.
   */ 

  myProps.get( dofTypes_, DOFS_PROP );

  if ( dofTypes_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "dofTypes must have the same length as nodeGroups" );
  }

  myConf.set ( DOFS_PROP, dofTypes_ );

  /*
   *  Get factors associated with the node groups. The length 
   *  of factors must be equal to the length of node groups.
   */ 

  myProps.get( factors_, FACTORS_PROP );

  if ( factors_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "factors must have the same length as nodeGroups" );
  }

  myConf.set ( FACTORS_PROP,   factors_      );

  // Get step size 

  stepSize_ = 0.0;

  myProps.find( stepSize_, STEP_SIZE_PROP );
  myConf. set ( STEP_SIZE_PROP, stepSize_ );

  // Set total and total0 to zero

  total_ = total0_ = 0.0;
  
}


PeriodicModel::~PeriodicModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void PeriodicModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  if ( ! props.contains( myName_ ) )
  {
    return;
  }

}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void PeriodicModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool PeriodicModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jem::Exception;
  using jive::model::Actions;
  using jive::model::ActionParams;

  // initialization

  if ( action == Actions::INIT )
  {
    init_ ( globdat );

    return true;
  }

  // apply displacement increment

  if ( action == Actions::GET_CONSTRAINTS )
  {
    applyConstraints_ ( params, globdat );

    return true;
  }

  // proceed to next time step

  if ( action == Actions::COMMIT )
  {
    commit_ ( params, globdat );

    return true;
  }

  // advance to next time step

  if ( action == Actions::ADVANCE )
  {
    globdat.set ( "var.accepted", true );

    advance_ ( globdat );

    return true;
  }

  // adapt step size

  else if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  return false;

}


//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void PeriodicModel::init_ ( const Properties& globdat )
{
  // Get nodes, then dofs of nodes, and constraints of dofs

  nodes_ = NodeSet::find    ( globdat );
  dofs_  = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat );

  idx_t                 nn0;        // number of nodes in group 0
  idx_t                 nn1;        // number of nodes in group 1
  Assignable<NodeGroup> group0;     
  Assignable<NodeGroup> group1;     
  IdxVector             inodes0;    // nodes in group 0
  IdxVector             inodes1;    // nodes in group 1
  String                context;

  // loop over node groups

  for ( idx_t ig = 0; ig < ngroups_; ig++ )
  {
    group0  = NodeGroup::get ( nodeGroups0_[ig], nodes_, globdat, context );
    group1  = NodeGroup::get ( nodeGroups1_[ig], nodes_, globdat, context );

    nn0     = group0.size();
    inodes0 . resize ( nn0 );
    inodes0 = group0.getIndices ();

    nn1     = group1.size();
    inodes1 . resize ( nn1 );
    inodes1 = group1.getIndices ();

    if ( nn0 != nn1 )
    {
      throw IllegalInputException ( JEM_FUNC,
          "periodic nodeGroups must have the same number of nodes" );
    }

    idx_t itype  = dofs_->findType ( dofTypes_[ig] );

    double rval = total_ * factors_[ig];
    double coef = 1.0;

    // apply constraint

    for ( idx_t in = 0; in < nn0; in++ )
    {
      idx_t jdof = dofs_->getDofIndex ( inodes0[in], itype );
      idx_t idof = dofs_->getDofIndex ( inodes1[in], itype );

      cons_->addConstraint ( idof, rval, jdof, coef );
    }
  }

  // compress for more efficient storage

  cons_->compress();

}


//-----------------------------------------------------------------------
//   advance_
//-----------------------------------------------------------------------


void PeriodicModel::advance_ ( const Properties& globdat )
{
  // Update the total value of Periodic rhs

  total_   = total0_ + stepSize_ ;

}


//-----------------------------------------------------------------------
//   applyConstraints_
//-----------------------------------------------------------------------

void PeriodicModel::applyConstraints_

  ( const Properties&  params,
    const Properties&  globdat )

{

  if ( abs(stepSize_) > 1.e-20 && sum ( factors_ ) > 1.e-20 )
  {

  idx_t                 nn1;        // number of nodes in group 1
  Assignable<NodeGroup> group1;     
  IdxVector             inodes1;    // nodes in group 1
  String                context;

  // loop over node groups

  for ( idx_t ig = 0; ig < ngroups_; ig++ )
  {
    group1  = NodeGroup::get ( nodeGroups1_[ig], nodes_, globdat, context );

    nn1     = group1.size();
    inodes1 . resize ( nn1 );
    inodes1 = group1.getIndices ();

    idx_t itype  = dofs_->findType ( dofTypes_[ig] );

    double rval = total_ * factors_[ig];

    // apply constraint

    for ( idx_t in = 0; in < nn1; in++ )
    {
      idx_t idof = dofs_->getDofIndex ( inodes1[in], itype );

      cons_->setRvalue ( idof, rval );
    }
  }

  // compress for more efficient storage

  cons_->compress();

  }

}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void PeriodicModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // Update total step

  total0_  = total_;
}


//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void PeriodicModel::setStepSize_

  ( const Properties&  params )

{
  double       dt;

  params.get ( dt,  XProps::STEP_SIZE   );

  System::out() << "Setting step size to " << dt << " ("
    << dt / stepSize_ * 100 << "\% of previous step size)" << "\n";

  stepSize_ = dt;
  
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newPeriodicModel
//-----------------------------------------------------------------------


Ref<Model>            newPeriodicModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<PeriodicModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declarePeriodicModel
//-----------------------------------------------------------------------


void declarePeriodicModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( PeriodicModel::TYPE_NAME, newPeriodicModel );
}
