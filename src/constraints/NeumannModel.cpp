
/** @file NeumannModel.cpp
 *  @brief Implements Neumann boundary conditions.
 *  
 *  This class implements the computation of traction
 *  loads on the Neumann boundary. It is compatible 
 *  with any FE model, since does not require any pre-
 *  knowledge of the dofs. The dofs are extracted from
 *  the (X)Dofspace. This model must be included after
 *  the FEmodel in the *.pro file.
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 13 March 2022
 *
 *  Usage:      
 * 
 *    model =
 *       {
 *        type     = "Neumann";
 *        elements = "TopElems";
 *        loads    = [0.0,-1.0e+04,0.0]; // size = dofs
 *        shape  =
 *         {
 *           type  = "BLine2";
 *           shapeFuncs  =
 *           {
 *             type  = "Linear";
 *           };
 *           intScheme = "Gauss1";
 *         };
 *       };
 *
 *  Updates (when, what and who)
 *     - [22 March 2022] throws an IllegalInputException
 *       when size of loads is not equal to dofs.
 *       
 */

/* Include jem and jive headers */

#include <jem/base/System.h>
#include <jem/base/Float.h>
#include <jem/numeric/utilities.h>
#include <jem/numeric/algebra/utilities.h>
#include <jive/model/Actions.h>

/* Include other headers */

#include "NeumannModel.h"

using jive::model::Actions;


//=======================================================================
//   class NeumannModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  NeumannModel::TYPE_NAME          = "Neumann";
const char*  NeumannModel::LOADS_PROP         = "loads";
const char*  NeumannModel::SHAPE_PROP         = "shape";


//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------

NeumannModel::NeumannModel

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  globdat ) :

    Super ( name )

{
  using jive::util::joinNames;
  using jive::geom::BShapeFactory;

  Properties    myProps = props.findProps ( myName_ );
  Properties    myConf  = props.makeProps ( myName_ );

  const String  context = getContext ();

  // Get element group assigned to this model.

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

  shape_ = BShapeFactory::newInstance (
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

  // Get the DOFSpace

  dofs_  = XDofSpace::get ( nodes_.getData(), globdat );

  // Get number of dofs from the XDofSpace

  const int          ndofs    = dofs_->typeCount   ();
  const StringVector dofNames = dofs_->getTypeNames();

  dofNames_.resize ( ndofs );
  dofTypes_.resize ( ndofs );
  
  dofNames_ = dofNames;

  for ( idx_t i = 0; i < ndofs; i++ )
  {
    dofTypes_[i] = dofs_->addType ( dofNames_[i] );
  }

  loads_   .resize ( ndofs );
  loads_   = 0.0;

  myProps.get ( loads_, LOADS_PROP );
  myConf. set ( LOADS_PROP, loads_ );

  if ( loads_.size() != ndofs )
  {
    throw IllegalInputException (
      context,
      String::format (
        "loads has invalid size: %d (should be %d)",
        loads_.size(),
        ndofs
      )
    );
  }
}


NeumannModel::~NeumannModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void NeumannModel::configure

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


void NeumannModel::getConfig

  ( const Properties&  conf,
    const Properties&  globdat ) const

{
  Properties  myConf = conf.makeProps ( myName_ );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool NeumannModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jem::Exception;
  using jive::model::Actions;
  using jive::model::ActionParams;

  if ( action == Actions::GET_EXT_VECTOR )
  {
    Vector fext;
    double scale;

    scale = 1.0;

    params.find ( scale, ActionParams::SCALE_FACTOR );
    params.get  ( fext,  ActionParams::EXT_VECTOR   );

    getExtForce_    ( fext, scale, globdat );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getExtForce_
//-----------------------------------------------------------------------

void NeumannModel::getExtForce_

  ( const Vector&       fext,
    double              scale,
    const Properties&   globdat ) const

  {
    using jem::select;

    IdxVector ielems    = egroup_.getIndices ();
    Matrix    sfuncs    = shape_->getShapeFunctions ();

    const int elemCount = ielems.size         ();
    const int nodeCount = shape_->nodeCount   ();
    const int ipCount   = shape_->ipointCount ();
    const int loadCount = dofTypes_.size      ();
    const int dofCount  = nodeCount * loadCount ;

    Vector    fe        ( dofCount          );
    IdxVector idofs     ( dofCount          );
    IdxVector inodes    ( nodeCount         );
    Vector    weights   ( ipCount           );
    Matrix    coords    ( rank_, nodeCount  );
    Matrix    normals   ( rank_, ipCount    );

    // Loop over surface elements

    for ( idx_t ie = 0; ie < elemCount; ie++ )
    {
      
      // Get the global element index.

      int ielem = ielems[ie];

      // Get the element coordinates and DOFs.

      elems_. getElemNodes  ( inodes, ielem  );
      nodes_. getSomeCoords ( coords, inodes );
      dofs_-> getDofIndices ( idofs,  inodes, dofTypes_ );
      shape_->getNormals    ( normals, weights, coords );

      // Initialize the element load vector

      fe = 0.0;

      // Loop on integration points 
    
      for ( int ip = 0; ip < ipCount; ip++ )
      {

        // Loop over the load dofs

        for ( int ldof = 0; ldof < loadCount; ldof++ )
        {
          double load = scale * weights[ip] * loads_[ldof];

          fe[slice(ldof,END,loadCount)] += load * sfuncs(ALL,ip);
        }
      }

      // assemble into the global external force vector

      select ( fext, idofs ) += fe;

    }
  }


//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newNeumannModel
//-----------------------------------------------------------------------


Ref<Model>            newNeumannModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<NeumannModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareNeumannModel
//-----------------------------------------------------------------------


void declareNeumannModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( NeumannModel::TYPE_NAME, newNeumannModel );
}
