/** @file DirichletModel.cpp
 *  @brief Implements Dirichlet boundary conditions.
 *  
 *  This class implements a model for setting Dirichlet 
 *  boundary conditions for certain group of nodes. It
 *  must be included after the FE Model in *.pro file.
 *  The step size set in this model may be modified by
 *  any module that throws the action 
 *  SolverNames::SET_STEP_SIZE.
 *
 *  Author: F.P. van der Meer, f.p.vandermeer@tudelft.nl
 *          R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: May 2010
 *  
 *  Usage:      
 * 
 *    model =
 *       {
 *        type       = "Dirichlet";
 *        nodeGroups = "TopNodes";
 *        dofs       = "dy";
 *        factors    = [ -1.];
 *        stepSize   = 0.001;
 *       };
 *
 *  Updates (when, what and who)
 *     - [28 March 2022] 
 */

/* Include jem and jive headers */

#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>

/* Include user-defined class headers */

#include "DirichletModel.h"
#include "util/XNames.h"

using jive::model::Actions;


//=======================================================================
//   class DirichletModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  DirichletModel::TYPE_NAME          = "Dirichlet";
const char*  DirichletModel::NODE_GROUP_PROP    = "nodeGroups";
const char*  DirichletModel::DOFS_PROP          = "dofs";
const char*  DirichletModel::FACTORS_PROP       = "factors";
const char*  DirichletModel::STEP_SIZE_PROP     = "stepSize";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


DirichletModel::DirichletModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  ngroups_  = 0;
  stepSize_ = 0.0;
  total_    = total0_ = 0.0;
}


DirichletModel::~DirichletModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void DirichletModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties    myProps = props.findProps ( myName_ );

  const String  context = getContext ();

  // Get the node groups assigned to this model

  myProps.get( nodeGroups_, NODE_GROUP_PROP );

  ngroups_ = nodeGroups_.size();

  // Get dofs associated with the node groups. The length of
  // dofs must be equal to the length of node groups. 

  myProps.get( dofTypes_, DOFS_PROP );

  if ( dofTypes_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "dofTypes must have the same length as nodeGroups" );
  }

  // Get factors associated with the node groups. The length 
  // of factors must be equal to the length of node groups.

  myProps.get( factors_, FACTORS_PROP );

  if ( factors_.size() != ngroups_ )
  {
    throw IllegalInputException ( JEM_FUNC,
          "factors must have the same length as nodeGroups" );
  }

  // Get step size 
  
  myProps.find( stepSize_, STEP_SIZE_PROP );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void DirichletModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );

  const String  context = getContext ();

  myConf.set ( NODE_GROUP_PROP, nodeGroups_ );
  myConf.set ( DOFS_PROP, dofTypes_         );
  myConf.set ( FACTORS_PROP,   factors_     );
  myConf.set ( STEP_SIZE_PROP, stepSize_    );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool DirichletModel::takeAction

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

  // apply constraints

  if ( action == Actions::GET_CONSTRAINTS )
  {
    applyConstraints_ ( params, globdat );

    return true;
  }

  // actions after accepting/committing to a solution

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

  // set step size

  else if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  return false;

}


//-----------------------------------------------------------------------
//   init_
//-----------------------------------------------------------------------


void DirichletModel::init_ ( const Properties& globdat )
{
  // Get nodes, then dofs of nodes, and constraints of dofs

  nodes_ = NodeSet::find    ( globdat );
  dofs_  = XDofSpace::get   ( nodes_.getData(), globdat );
  cons_  = Constraints::get ( dofs_, globdat );    
}


//-----------------------------------------------------------------------
//   advance_
//-----------------------------------------------------------------------


void DirichletModel::advance_ ( const Properties& globdat )
{
  // Update the total value of Dirichlet rhs

  total_   = total0_ + stepSize_ ;

}


//-----------------------------------------------------------------------
//   applyConstraints_
//-----------------------------------------------------------------------

void DirichletModel::applyConstraints_

  ( const Properties&  params,
    const Properties&  globdat )

{
  idx_t                 nn;
  Assignable<NodeGroup> group;
  IdxVector             inodes;
  String                context;

  // loop over node groups

  for ( idx_t ig = 0; ig < ngroups_; ig++ )
  {
    group  = NodeGroup::get ( nodeGroups_[ig], nodes_, globdat, context );

    nn     = group.size();
    inodes . resize ( nn );
    inodes = group.getIndices ();

    idx_t itype  = dofs_->findType ( dofTypes_[ig] );

    double val = total_ * factors_[ig];

    // apply constraint

    for ( idx_t in = 0; in < nn; in++ )
    {
      idx_t idof = dofs_->getDofIndex ( inodes[in], itype );

      cons_->addConstraint ( idof, val );
    }
  }

  // compress for more efficient storage

  cons_->compress();
}


//-----------------------------------------------------------------------
//   commit_
//-----------------------------------------------------------------------


void DirichletModel::commit_

  ( const Properties&  params,
    const Properties&  globdat )

{
  // Update total step

  total0_  = total_;
}


//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void DirichletModel::setStepSize_

  ( const Properties&  params )

{
  double       dt;

  params.get ( dt, XProps::STEP_SIZE );

  System::out() << "Setting step size to " << dt << " ("
    << dt / stepSize_ * 100 << "\% of previous step size)" << "\n";

  stepSize_ = dt;
  
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newDirichletModel
//-----------------------------------------------------------------------


Ref<Model>            newDirichletModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<DirichletModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareDirichletModel
//-----------------------------------------------------------------------


void declareDirichletModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( DirichletModel::TYPE_NAME, newDirichletModel );
}
