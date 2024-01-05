
/** @file PoissonModel.cpp
 *  @brief Implements a Poisson model.
 *
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 11 November 2022
 *
 *  Updates (when, what and who)
 *     - [XX YYYY 2022] 
 * 
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>

/* Include other headers */

#include "FalconBasicModels.h"
#include "PoissonModel.h"
#include "util/TbFiller.h"


//=======================================================================
//   class PoissonModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  PoissonModel::DOF_NAME[1]      = { "u" };
const char*  PoissonModel::SHAPE_PROP       = "shape";
const char*  PoissonModel::COEFF_PROP       = "a";
const char*  PoissonModel::SOURCE_PROP      = "source";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


PoissonModel::PoissonModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::geom::IShapeFactory;

  Properties    myProps = props.getProps  ( myName_ );
  Properties    myConf  = conf .makeProps ( myName_ );

  const String  context = getContext ();

  // Get the element group assigned to this model.

  egroup_ = ElemGroup::get ( myConf, myProps, globdat, context );

  elems_  = egroup_.getElements ();
  nodes_  = elems_ .getNodes    ();
  rank_   = nodes_ .rank        ();

  // Make sure that the number of spatial dimensions (the rank of the
  // mesh) is valid.

  if ( rank_ < 1 || rank_ > 3 )
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  }

  // Create an internal shape object for computing the element shape
  // functions.

  shape_ = IShapeFactory::newInstance (
    joinNames ( myName_, SHAPE_PROP ),
    conf,
    props
  );

  // Make sure that the rank of the shape matches the rank of the
  // mesh.

  if ( shape_->globalRank() != rank_ )
  {
    throw IllegalInputException (
      context,
      String::format (
        "shape has invalid rank: %d (should be %d)",
        shape_->globalRank (),
        rank_
      )
    );
  }

  // Make sure that each element has the same number of nodes as the
  // shape object.

  elems_.checkSomeElements (
    context,
    egroup_.getIndices (),
    shape_->nodeCount  ()
  );


  // Get the DOF space, add displacement and pressure DOFs

  dofs_ = XDofSpace::get ( nodes_.getData(), globdat );

  dofTypes_.resize ( 1 );
  dofTypes_[0] = dofs_->addType ( DOF_NAME[0] );
  
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  a_       = 1.0;
  f_       = 1.0;
 
}


PoissonModel::~PoissonModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void PoissonModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props  .findProps ( myName_ );

  myProps.find ( a_, COEFF_PROP  );
  myProps.find ( f_, SOURCE_PROP );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void PoissonModel::getConfig ( const Properties& conf,
                               const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );

  myConf. set  ( COEFF_PROP,  a_ );
  myConf. set  ( SOURCE_PROP, f_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool PoissonModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {

    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  force;

    // Get the current state.

    StateVector::get    ( state, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0 );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, state );

    return true;
  }

  // things to be done when the solution converges, ie the end of 
  // a load step

  if ( action == Actions::COMMIT )
  {
    return true;
  }

  // used with the StepModule
  // to advance to new load step or resolve with the same load vector

  if ( action == Actions::CHECK_COMMIT )
  {
    return true;
  }

  // post processing: compute some table for visualize stress, 
  // strain, ...

  if ( action == Actions::GET_TABLE )
  {
    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void PoissonModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   state )

{

  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state
  
  Vector      elemState  ( nodeCount );

  // internal force vector

  Vector      elemForce  ( nodeCount );

  // element stiffness matrix

  Matrix      elemMat    ( nodeCount, nodeCount );

  // B matrices

  Matrix      be         ( rank_, nodeCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   idofs      ( nodeCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( idofs,    inodes, dofTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current nodal state

    elemState = select ( state, idofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce  = 0.0;
    elemMat    = 0.0;    

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // get B-matrix
      
      be =  grads(ALL,ALL,ip);

      Vector Nip          = N( ALL,ip );

      // compute the stiffness matrix

      wip         = ipWeights[ip];
      elemMat    += wip * a_ * mc2.matmul ( bet, be );
     
      // compute internal forces

      elemForce  +=  wip * ( a_ * mc1.matmul( mc2.matmul ( bet, be ), elemState ) - Nip * f_ );

    }  // end of loop on integration points

    // Assembly ...

    mbuilder.addBlock ( idofs, idofs, elemMat );

    select ( force, idofs )  += elemForce;
  }
}


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newPoissonModel
//-----------------------------------------------------------------------


static Ref<Model>     newPoissonModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<PoissonModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declarePoissonModel
//-----------------------------------------------------------------------


void declarePoissonModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "Poisson", & newPoissonModel );
}
