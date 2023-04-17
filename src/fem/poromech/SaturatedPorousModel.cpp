/** @file SaturatedPorousModel.cpp
 *  @brief Implements saturated porous media model.
 *  
 *  This class implements the fully saturated porous media 
 *  model with two phases - solid and fluid. The effect of
 *  gravity is absent. The fluid pressure equation is 
 *  stabilized via perturbation (see DOI: 10.1002/nme2295).
 *  
 *  Author: R. Bharali, ritukesh.bharali@chalmers.se
 *  Date: 04 March 2022
 * 
 *  NOTE: Stiffness matrix is unsymmetric, choose the linear
 *        solver accordingly. Stabilization works only for 
 *        Triangle3 elements.
 * 
 *  TO-DO: Anisotropic permeability
 *
 *  Updates (when, what and who)
 *     - [11 March 2022] added option for perturbation 
 *       based numerical stabilization (RB), 
 *       DOI: 10.1002/nme2295
 * 
 *     - [14 March 2022] modified pressure equation to
 *       obtain a symmetric stiffness matrix (RB)
 * 
 *     - [22 March 2022] time step-size (dtime) can be
 *       modified by the solver that throws an action
 *       XActions::SET_STEP_SIZE. An example would
 *       be the AdaptiveSteppingModule. (RB)
 * 
 *     - [05 July 2022] added getIntForce_ to compute
 *       internal force, required by quasi-Newton/line
 *       search algorithms. (RB)
 *       
 */

/* Include jem and jive headers */

#include <jem/base/Error.h>
#include <jem/base/IllegalInputException.h>
#include <jive/model/Actions.h>
#include <jive/geom/Geometries.h>

/* Include other headers */

#include "FalconPoroMechModels.h"
#include "SaturatedPorousModel.h"
#include "util/TbFiller.h"
#include "util/XNames.h"

#include "materials/Material.h"

//=======================================================================
//   class SaturatedPorousModel
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  SaturatedPorousModel::DISP_NAMES[3]     = { "dx", "dy", "dz" };
const char*  SaturatedPorousModel::SHAPE_PROP        = "shape";
const char*  SaturatedPorousModel::MATERIAL_PROP     = "material";

const char*  SaturatedPorousModel::INTRIN_PERM_PROP  = "intrin_perm";
const char*  SaturatedPorousModel::FLUID_VISC_PROP   = "fluid_visc";
const char*  SaturatedPorousModel::SOLID_STIFF_PROP  = "solid_stiff";
const char*  SaturatedPorousModel::FLUID_STIFF_PROP  = "fluid_stiff";
const char*  SaturatedPorousModel::POROSITY_PROP     = "porosity";
const char*  SaturatedPorousModel::BIOT_COEFF_PROP   = "biot_coeff";
const char*  SaturatedPorousModel::DTIME_PROP        = "dtime";

//-----------------------------------------------------------------------
//   constructors & destructor
//-----------------------------------------------------------------------


SaturatedPorousModel::SaturatedPorousModel

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

  dofTypes_.resize ( rank_ + 1 );

  // Make a shallow copy of the first part of the dofTypes_ array.

  dispTypes_.ref ( dofTypes_[slice(BEGIN,rank_)] );

  for ( int i = 0; i < rank_; i++ )
  {
    dispTypes_[i] = dofs_->addType ( DISP_NAMES[i] );
  }

  // Make a shallow copy of the last part of the dofTypes_ array.

  presTypes_.ref ( dofTypes_[slice(rank_,END)] );

  presTypes_[0] = dofs_->addType ( "dp" );
    
  dofs_->addDofs (
    elems_.getUniqueNodesOf ( egroup_.getIndices() ),
    dofTypes_
  );

  // Compute the total number of integration points.

  int  ipCount = shape_->ipointCount() * egroup_.size();

  // Allocate memory for the strains and set them to zero

  strain_.resize ( STRAIN_COUNTS[rank_], ipCount );
  strain_ = 0.0;

  // Create a material model object.

  material_ = newMaterial ( MATERIAL_PROP, myConf, myProps, globdat );

  material_->allocPoints  ( ipCount );

  initializeIPMPMap_ ();

  // Select the correct function for computing the B-matrix.

  getShapeGrads_ = getShapeGradsFunc ( rank_ );

  IdxVector   ielems     = egroup_.getIndices ();
  const int   ielemCount = ielems.size         ();

  isActive_.resize ( ielemCount );
  isActive_ = 1;

  // Configure poromechanical parameters

  kappa_       = 1.0e-12;
  myProps.find ( kappa_, INTRIN_PERM_PROP );
  myConf. set  ( INTRIN_PERM_PROP, kappa_ );

  mu_          = 8.9e-4;
  myProps.find ( mu_, FLUID_VISC_PROP );
  myConf. set  ( FLUID_VISC_PROP, mu_ );

  Ks_          = 1.0e+10;
  myProps.find ( Ks_, SOLID_STIFF_PROP );
  myConf. set  ( SOLID_STIFF_PROP, Ks_ );

  Kf_          = 2.0e+09;
  myProps.find ( Kf_, FLUID_STIFF_PROP );
  myConf. set  ( FLUID_STIFF_PROP, Kf_ );

  phi_         = 0.375;
  myProps.find ( phi_, POROSITY_PROP );
  myConf. set  ( POROSITY_PROP, phi_ );

  alpha_      = 1.0;
  myProps.find ( alpha_, BIOT_COEFF_PROP );
  myConf. set  ( BIOT_COEFF_PROP, alpha_ );

  dtime_       = 1.0;
  myProps.find ( dtime_, DTIME_PROP );
  myConf. set  ( DTIME_PROP, dtime_ );

  // Compute derived quantities

  Sto_  = ( phi_ / Kf_ ) + ( alpha_ - phi_ ) / Ks_;
  Keff_ = kappa_ / mu_ ;

  // Set up Voigt divergence operator based of mesh rank

  voigtDiv_.resize ( STRAIN_COUNTS[rank_]);
  voigtDiv_ = 0.0;

  if ( rank_ == 1 )
  {
    voigtDiv_[0] = 1.0;      // xx-component
  }
  else if ( rank_ == 2 )
  {
    voigtDiv_[0] = 1.0;      // xx-component
    voigtDiv_[1] = 1.0;      // yy-component
  }
  else if ( rank_ == 3 )
  {
    voigtDiv_[0] = 1.0;      // xx-component
    voigtDiv_[1] = 1.0;      // yy-component
    voigtDiv_[2] = 1.0;      // zz-component
  }
  else
  {
    throw IllegalInputException (
      context,
      String::format (
        "invalid node rank: %d (should be 1, 2 or 3)", rank_
      )
    );
  } 
}


SaturatedPorousModel::~SaturatedPorousModel ()
{}


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void SaturatedPorousModel::configure

  ( const Properties&  props,
    const Properties&  globdat )

{
  Properties  myProps  = props  .findProps ( myName_ );
  Properties  matProps = myProps.findProps ( MATERIAL_PROP );

  material_->configure ( matProps, globdat );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void SaturatedPorousModel::getConfig ( const Properties& conf,
                                      const Properties&  globdat )  const
{
  Properties  myConf  = conf  .makeProps ( myName_ );
  Properties  matConf = myConf.makeProps ( MATERIAL_PROP );

  material_->getConfig ( matConf, globdat );
}


//-----------------------------------------------------------------------
//   takeAction
//-----------------------------------------------------------------------


bool SaturatedPorousModel::takeAction

  ( const String&      action,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  // compute the internal force vector

  if ( action == Actions::GET_INT_VECTOR )
  {

    Vector  state;
    Vector  state0;
    Vector  force;

    // Get the current and old step displacements.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );

    // Get the internal force vector.

    params.get ( force,    ActionParams::INT_VECTOR );

    getIntForce_ ( force, state, state0 );

    return true;
  }

  // compute the tangent stiffness matrix and internal force vector

  if ( action == Actions::GET_MATRIX0 )
  {

    Ref<MatrixBuilder>  mbuilder;

    Vector  state;
    Vector  state0;
    Vector  force;

    // Get the current and old step displacements.

    StateVector::get    ( state,  dofs_, globdat );
    StateVector::getOld ( state0, dofs_, globdat );

    // Get the matrix builder and the internal force vector.

    params.get ( mbuilder, ActionParams::MATRIX0    );
    params.get ( force,    ActionParams::INT_VECTOR );

    getMatrix_ ( *mbuilder, force, state, state0 );

    return true;
  }

  // compute mass matrix 

  if ( action == Actions::GET_MATRIX2 )
  {
    Ref<MatrixBuilder> mbuilder;

    params.get ( mbuilder, ActionParams::MATRIX2 );
    
    getMatrix2_( *mbuilder );

    return true;
  }

  // During initialization of the model, the pressure dofs of the
  // mid nodes of all domain elements are linearly constrained to
  // the corner nodes of the respective edge. This yields a Taylor
  // Hood element. 

  if ( action == Actions::INIT )
  {
    // Get the constraints associated with the DOF space.

    Ref<Constraints>  cons = Constraints::get ( dofs_, globdat );

    setConstraints_ ( *cons );

    return true;
  }

  // COMMIT: Requests an action when a computation step has converged.
  // Ideally, at this point, the material internal state variables
  // are swapped. 

  if ( action == Actions::COMMIT )
  {
    material_->commit ();

    return true;
  }

  // CHECK_COMMIT: Can used to discard the current step 

  if ( action == Actions::CHECK_COMMIT )
  {
    //checkCommit_ ( params );
    return true;
  }

  // GET_TABLE: Requests a post-processing action, such as computing
  // tabular data (stress, strains, etc.) to print VTK files. 
  
  if ( action == Actions::GET_TABLE )
  {
    return getTable_ ( params, globdat );
  }

  if ( action == XActions::SET_STEP_SIZE )
  {
    setStepSize_ ( params );
  }

  return false;
}


//-----------------------------------------------------------------------
//   getIntForce_
//-----------------------------------------------------------------------


void SaturatedPorousModel::getIntForce_

  ( const Vector&   force,
    const Vector&   state,
    const Vector&   state0 )

{
  // System::out() << "Enters getIntForce_ \n";

  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of displacement and pressure dofs

  const int   dispCount  = nodeCount * rank_;  
  const int   presCount   = nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement and pressure
  
  Vector      pres       ( presCount );
  Vector      disp       ( dispCount );

  // old step element vector state:
  // displacement and pressure
  
  Vector      pres0      ( presCount );
  Vector      disp0      ( dispCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );

  // old-step

  Vector      strain0    ( strCount );

  // internal force vector:
  // displacement and pressure

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( presCount  );

  // element stiffness matrices (four components)

  // Matrix      elemMat1   ( dispCount, dispCount ); // disp-disp
  Matrix      elemMat2   ( dispCount, presCount ); // disp-pres
  Matrix      elemMat3   ( presCount, dispCount ); // pres-disp
  Matrix      elemMat4   ( presCount, presCount ); // pres-pres
  
  // B matrices
  // displacement and pressure
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, presCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   presDofs   ( presCount );
  IdxVector   dispDofs   ( dispCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {

    if ( ! isActive_[ie] )
    {      
      continue;
    }

    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( presDofs, inodes, presTypes_ );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current nodal displacements and pressures

    pres = select ( state, presDofs );
    disp = select ( state, dispDofs );

    // Get old step nodal displacements and pressures

    pres0 = select ( state0, presDofs );
    disp0 = select ( state0, dispDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0;

    // elemMat1   = 0.0;
    elemMat2   = 0.0;
    elemMat3   = 0.0;
    elemMat4   = 0.0;    

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {

      // get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );

      // Compute the B-matrix for this integration point.
      // it is the B-matrix of displacement dofs

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // get B-matrix associated with pres dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the strain, strain0 for this integration point.

      matmul ( strain,  bd, disp  );
      matmul ( strain0, bd, disp0 );

      // Store the regular strain components.

      strain_(ALL,ipoint) = strain;

      // Compute the integration point pressures

      Vector Nip          = N( ALL,ip );

      double p            = dot ( Nip, pres  );
      double p0           = dot ( Nip, pres0 );

      // Update the material model.

      material_->update ( stress, stiff, strain_(ALL,ipoint), ipoint );

      // compute stiffness matrix components

      wip         = ipWeights[ip];
      // elemMat1   += wip * mc3.matmul ( bdt, stiff, bd );
      elemMat2   -= wip * alpha_ * mc2.matmul ( bdt, matmul (voigtDiv_, Nip ) );
      elemMat4   -= wip * ( Sto_ * matmul ( Nip, Nip ) + dtime_ * Keff_ * mc2.matmul ( bet, be ) );

      /** (Non-symmetric) 
       * elemMat4   += wip * ( Sto_ * matmul ( Nip, Nip ) 
       *  + dtime_  * Keff_ * mc2.matmul ( bet, be ) ); */
     
      // compute internal forces

      elemForce1 +=  wip * ( mc1.matmul ( bdt, stress ) );
      elemForce1 -=  wip * alpha_ * p * ( mc1.matmul ( bdt, voigtDiv_ ) );
      elemForce2 +=  wip * ( Sto_ * p0 * Nip );

      /** (Non-symmetric) 
       * elemForce2 -=  wip * ( Sto_ * p0 * Nip ); */

    }  // end of loop on integration points

    /** (Non-symmetric) 
       * elemMat3      =  -elemMat2.transpose() ; */

    elemMat3      =  elemMat2.transpose() ;
    elemForce2   +=  mc1.matmul ( elemMat3, disp ) - mc1.matmul ( elemMat3, disp0 );
    elemForce2   +=  mc1.matmul ( elemMat4, pres );

    // Assembly ...

    select ( force, dispDofs )  += elemForce1;
    select ( force, presDofs  ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getMatrix_
//-----------------------------------------------------------------------


void SaturatedPorousModel::getMatrix_

  ( MatrixBuilder&  mbuilder,
    const Vector&   force,
    const Vector&   state,
    const Vector&   state0 )

{
  // System::out() << "Enters getMatrix_ \n";

  IdxVector   ielems     = egroup_.getIndices ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   strCount   = STRAIN_COUNTS[rank_];

  // number of displacement and pressure dofs

  const int   dispCount  = nodeCount * rank_;  
  const int   presCount   = nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );

  // current element vector state:
  // displacement and pressure
  
  Vector      pres       ( presCount );
  Vector      disp       ( dispCount );

  // old step element vector state:
  // displacement and pressure
  
  Vector      pres0      ( presCount );
  Vector      disp0      ( dispCount );

  // current

  Matrix      stiff      ( strCount, strCount );
  Vector      stress     ( strCount );
  Vector      strain     ( strCount );

  // old-step

  Vector      strain0    ( strCount );

  // internal force vector:
  // displacement and pressure

  Vector      elemForce1 ( dispCount );
  Vector      elemForce2 ( presCount  );

  // element stiffness matrices (four components)

  Matrix      elemMat1   ( dispCount, dispCount ); // disp-disp
  Matrix      elemMat2   ( dispCount, presCount ); // disp-pres
  Matrix      elemMat3   ( presCount, dispCount ); // pres-disp
  Matrix      elemMat4   ( presCount, presCount ); // pres-pres
  
  // B matrices
  // displacement and pressure
  
  Matrix      bd         ( strCount, dispCount );
  Matrix      bdt        = bd.transpose ();

  Matrix      be         ( rank_, presCount );
  Matrix      bet        = be.transpose ();
  
  IdxVector   inodes     ( nodeCount );
  IdxVector   presDofs   ( presCount );
  IdxVector   dispDofs   ( dispCount );

  Vector      ipWeights  ( ipCount   );

  MChain1     mc1;
  MChain2     mc2;
  MChain3     mc3;
  MChain4     mc4; 
  
  double      wip;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {

    if ( ! isActive_[ie] )
    {      
      continue;
    }

    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes,   ielem  );
    nodes_.getSomeCoords ( coords,   inodes );
    dofs_->getDofIndices ( presDofs, inodes, presTypes_ );
    dofs_->getDofIndices ( dispDofs, inodes, dispTypes_ );

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Compute the shape functions and transpose

    Matrix N  = shape_->getShapeFunctions ();
    
    // Get current nodal displacements and pressures

    pres = select ( state, presDofs );
    disp = select ( state, dispDofs );

    // Get old step nodal displacements and pressures

    pres0 = select ( state0, presDofs );
    disp0 = select ( state0, dispDofs );

    // Initialize the internal forces
    // and the tangent stiffness matrix
    
    elemForce1 = 0.0;
    elemForce2 = 0.0;

    elemMat1   = 0.0;
    elemMat2   = 0.0;
    elemMat3   = 0.0;
    elemMat4   = 0.0;    

    // Loop on integration points 
    
    for ( int ip = 0; ip < ipCount; ip++ )
    {

      // get the corresponding material point ipoint
      // of element ielem integration point ip from the mapping

      int ipoint = ipMpMap_ ( ielem, ip );

      // Compute the B-matrix for this integration point.
      // it is the B-matrix of displacement dofs

      getShapeGrads_ ( bd, grads(ALL,ALL,ip) );

      // get B-matrix associated with pres dofs
      
      be =  grads(ALL,ALL,ip);

      // Compute the strain, strain0 for this integration point.

      matmul ( strain,  bd, disp  );
      matmul ( strain0, bd, disp0 );

      // Store the regular strain components.

      strain_(ALL,ipoint) = strain;

      // Compute the integration point pressures

      Vector Nip          = N( ALL,ip );

      double p            = dot ( Nip, pres  );
      double p0           = dot ( Nip, pres0 );

      // Update the material model.

      material_->update ( stress, stiff, strain_(ALL,ipoint), ipoint );

      // compute stiffness matrix components

      wip         = ipWeights[ip];
      elemMat1   += wip * mc3.matmul ( bdt, stiff, bd );
      elemMat2   -= wip * alpha_ * mc2.matmul ( bdt, matmul (voigtDiv_, Nip ) );
      elemMat4   -= wip * ( Sto_ * matmul ( Nip, Nip ) + dtime_ * Keff_ * mc2.matmul ( bet, be ) );

      /** (Non-symmetric) 
       * elemMat4   += wip * ( Sto_ * matmul ( Nip, Nip ) 
       *  + dtime_  * Keff_ * mc2.matmul ( bet, be ) ); */
     
      // compute internal forces

      elemForce1 +=  wip * ( mc1.matmul ( bdt, stress ) );
      elemForce1 -=  wip * alpha_ * p * ( mc1.matmul ( bdt, voigtDiv_ ) );
      elemForce2 +=  wip * ( Sto_ * p0 * Nip  );

      /** (Non-symmetric) 
       * elemForce2 -=  wip * ( Sto_ * p0 * Nip ); */

    }  // end of loop on integration points

    /** (Non-symmetric) 
       * elemMat3      =  -elemMat2.transpose() ; */

    elemMat3      =  elemMat2.transpose() ;
    elemForce2   +=  mc1.matmul ( elemMat3, disp ) - mc1.matmul ( elemMat3, disp0 );
    elemForce2   +=  mc1.matmul ( elemMat4, pres );

    // Assembly ...

    mbuilder.addBlock ( dispDofs, dispDofs, elemMat1 );
    mbuilder.addBlock ( dispDofs, presDofs, elemMat2 );
    mbuilder.addBlock ( presDofs, dispDofs, elemMat3 );
    mbuilder.addBlock ( presDofs, presDofs, elemMat4 );

    select ( force, dispDofs )  += elemForce1;
    select ( force, presDofs  ) += elemForce2;
  }
}


//-----------------------------------------------------------------------
//   getMatrix2_
//-----------------------------------------------------------------------

// compute the mass matrix 
// current implementation: consistent mass matrix

void SaturatedPorousModel::getMatrix2_

    ( MatrixBuilder&          mbuilder )
{
  IdxVector   ielems     = egroup_.getIndices  ();

  const int   ielemCount = ielems.size         ();
  const int   nodeCount  = shape_->nodeCount   ();
  const int   ipCount    = shape_->ipointCount ();
  const int   dofCount   = rank_ * nodeCount;

  Cubix       grads      ( rank_, nodeCount, ipCount );
  Matrix      coords     ( rank_, nodeCount );
  
  Matrix      elemMat    ( dofCount, dofCount );

  Matrix      R          ( rank_, rank_ );

  Matrix      sfuncs     = shape_->getShapeFunctions ();
  Matrix      N          ( rank_, rank_ * nodeCount );
  Matrix      Nt         = N.transpose ( ); 

  IdxVector   inodes     ( nodeCount );
  IdxVector   idofs      ( dofCount  );

  Vector      ipWeights  ( ipCount   );

  MChain3     mc3;

  double      rho = 0.0;

  R = 0.0;
 
  for ( int i = 0; i < rank_ ; i++ )
  {
    R(i,i) = rho;
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element coordinates and DOFs.

    elems_.getElemNodes  ( inodes, ielem    );
    nodes_.getSomeCoords ( coords, inodes );
    dofs_->getDofIndices ( idofs,  inodes, dispTypes_ );    // dofTypes_

    // Compute the spatial derivatives of the element shape
    // functions in the integration points.

    shape_->getShapeGradients ( grads, ipWeights, coords );

    // Assemble the element matrix and the internal force vector.

    elemMat   = 0.0;

    for ( int ip = 0; ip < ipCount; ip++ )
    {
      // compute matrix of shpae function N       

      getShapeFuncs_ ( N, sfuncs(ALL,ip) );

      // Add the contribution of this integration point.

      elemMat   += ipWeights[ip] * mc3.matmul ( Nt, R, N );
    }

    // Add the element secant matrix to the global stiffness matrix.

    mbuilder.addBlock ( idofs, idofs, elemMat );

  }
}

//-----------------------------------------------------------------------
//   setConstraints_
//-----------------------------------------------------------------------


void SaturatedPorousModel::setConstraints_ ( Constraints& cons )
{
  using jive::util::Printer;
  using jive::geom::Geometries;

  IdxVector     ielems     = egroup_.getIndices ();

  const int     ielemCount = ielems   .size     ();
  const int     nodeCount  = shape_->nodeCount  ();
  const String  shapeGeom  = shape_->getGeometry ();

  IdxVector     inodes     ( nodeCount );
  IdxVector     idofs      ( nodeCount );
  IdxVector     jdofs      ( 2 );
  Vector        coeffs     ( 2 );

  IntMatrix     cdofs;


  // only do this constraint for high order elements
  // not for linear elements

  if      ( shapeGeom == Geometries::LINE &&
      nodeCount == 2 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::TRIANGLE &&
      nodeCount == 3 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::SQUARE &&
      nodeCount == 4 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::TETRAHEDRON &&
      nodeCount == 4 )
  {
    return;
  }
  else if ( shapeGeom == Geometries::CUBE &&
      nodeCount == 8 )
  {
    return;
  }

  System::info() << " Initializing constraints for Taylor-Hood element ...\n";

  // Determine which pressure DOFS are to be constrained. The
  // first row of the cdofs matrix contains the DOF indices to be
  // constrained, and the other two rows contain the DOF indices of
  // the two master nodes. Note that the DOF indices are local with
  // respect to an element. They will be translated to global DOF
  // indices below.

  if      ( shapeGeom == Geometries::LINE &&
      nodeCount == 3 )
  {
    cdofs.resize ( 3, 1 );

    // Constrain the mid node.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;
  }
  else if ( shapeGeom == Geometries::TRIANGLE &&
      nodeCount == 6 )
  {
    cdofs.resize ( 3, 3 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 0;
  }
  else if ( shapeGeom == Geometries::SQUARE &&
      nodeCount == 8 )
  {
    cdofs.resize ( 3, 4 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 6;

    cdofs(0,3) = 7;
    cdofs(1,3) = 6;
    cdofs(2,3) = 0;
  }
  else if ( shapeGeom == Geometries::TETRAHEDRON &&
      nodeCount == 10 )
  {
    cdofs.resize ( 3, 6 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 0;

    cdofs(0,3) = 6;
    cdofs(1,3) = 9;
    cdofs(2,3) = 0;

    cdofs(0,4) = 7;
    cdofs(1,4) = 9;
    cdofs(2,4) = 1;

    cdofs(0,5) = 8;
    cdofs(1,5) = 9;
    cdofs(2,5) = 4;
  }
  else if ( shapeGeom == Geometries::CUBE &&
      nodeCount == 20 )
  {
    cdofs.resize ( 3, 12 );

    // Constrain the mid nodes on each edge.

    cdofs(0,0) = 1;
    cdofs(1,0) = 0;
    cdofs(2,0) = 2;

    cdofs(0,1) = 3;
    cdofs(1,1) = 2;
    cdofs(2,1) = 4;

    cdofs(0,2) = 5;
    cdofs(1,2) = 4;
    cdofs(2,2) = 6;

    cdofs(0,3) = 7;
    cdofs(1,3) = 6;
    cdofs(2,3) = 0;

    cdofs(0,4) = 8;
    cdofs(1,4) = 0;
    cdofs(2,4) = 12;

    cdofs(0,5) = 9;
    cdofs(1,5) = 2;
    cdofs(2,5) = 14;

    cdofs(0,6) = 10;
    cdofs(1,6) = 4;
    cdofs(2,6) = 16;

    cdofs(0,7) = 11;
    cdofs(1,7) = 6;
    cdofs(2,7) = 18;

    cdofs(0,8) = 15;
    cdofs(1,8) = 14;
    cdofs(2,8) = 16;

    cdofs(0,9) = 17;
    cdofs(1,9) = 16;
    cdofs(2,9) = 18;

    cdofs(0,10) = 19;
    cdofs(1,10) = 18;
    cdofs(2,10) = 12;

    cdofs(0,11) = 13;
    cdofs(1,11) = 12;
    cdofs(2,11) = 14;
  }
  else
  {
    throw IllegalInputException (
      getContext (),
      String::format (
  "unsupported shape geometry: %S%d",
  &shapeGeom,
  nodeCount
      )
    );
  }

  // Determine the constraint coefficients.

  coeffs = 0.5;

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element DOFs.

    elems_.getElemNodes  ( inodes, ielem );
    dofs_->getDofIndices ( idofs,  inodes, presTypes_ );

    // Add constraints.

    for ( int j = 0; j < cdofs.size(1); j++ )
    {
      // Get the global index of the DOF to be constrained.

      int  idof = idofs[cdofs(0,j)];

      // Check whether this DOF has already been constrained.

      if ( cons.isSlaveDof( idof ) )
      {
        //continue;
        cons.eraseConstraint( idof );
      }

      // Get the global indices of the two master DOFs.

      jdofs[0] = idofs[cdofs(1,j)];
      jdofs[1] = idofs[cdofs(2,j)];

      cons.addConstraint ( idof, jdofs, coeffs );

    }
  }
}

//-----------------------------------------------------------------------
//   initializeIPMPMap_
//-----------------------------------------------------------------------


void  SaturatedPorousModel::initializeIPMPMap_ ( )

{
  jive::IdxVector   ielems     = egroup_.getIndices  ();

  const idx_t   ielemCount = ielems.size         ();
  const idx_t   ipCount    = shape_->ipointCount ();

        idx_t   ipoint     = 0;

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // loop over integration points 

    for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
    {
      ipMpMap_ ( ielem, ip ) = ipoint;
    }
  }
}


//-----------------------------------------------------------------------
//   getTable_
//-----------------------------------------------------------------------


bool SaturatedPorousModel::getTable_

  ( const Properties&  params,
    const Properties&  globdat )

{
  using jive::model::Actions;
  using jive::model::ActionParams;
  using jive::model::StateVector;

  Ref<XTable>  table;
  Vector       weights;
  String       name;


  // Get the table, the name of the table, and the table row weights
  // from the action parameters.

  params.get ( table,   ActionParams::TABLE );
  params.get ( name,    ActionParams::TABLE_NAME );
  params.get ( weights, ActionParams::TABLE_WEIGHTS );

  // Nodal stresses
  if ( name == "stress" && 
       table->getRowItems() == nodes_.getData() )
  {
    Vector disp;

    StateVector::get (disp, dofs_, globdat);

    getStress_ ( *table, weights);

    return true;
  }

  // Nodal strains
  if ( name == "strain" && 
       table->getRowItems() == nodes_.getData() )
  {
    Vector disp;

    StateVector::get (disp, dofs_, globdat);

    getStrain_ ( *table, weights);

    return true;
  }

  // 
  if ( name == "xoutTable" )
  {
    Vector  disp;
    String  contents;

    StateVector::get ( disp,     dofs_, globdat  );
    params.      get ( contents, "contentString" );

    
    getOutputData_ ( table, weights, contents, disp );

    return true;
  }

  return false;
}


//-----------------------------------------------------------------------
//   getStress_
//-----------------------------------------------------------------------


void SaturatedPorousModel::getStress_

  ( XTable&        table,
    const Vector&  weights )

{
  Matrix     sfuncs     = shape_->getShapeFunctions ();

  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];
  // const int  dispCount   = nodeCount * rank_;

  

  Matrix     ndStress   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  Matrix     stiff      ( strCount, strCount );
  Vector     stressIp   ( strCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  MChain1    mc1;

  int        ipoint, igpoint = 0;

  // Add the columns for the stress components to the table.

  switch ( strCount )
  {
  case 1:

    jcols[0] = table.addColumn ( "stress_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_xy" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );

    break;  

  case 6:

    jcols[0] = table.addColumn ( "stress_xx" );
    jcols[1] = table.addColumn ( "stress_yy" );
    jcols[2] = table.addColumn ( "stress_zz" );
    jcols[3] = table.addColumn ( "stress_xy" );
    jcols[4] = table.addColumn ( "stress_yz" );
    jcols[5] = table.addColumn ( "stress_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of stress components: " +
      String ( strCount )
    );
  }

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndStress  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount; ip++, igpoint++ )
    {
      // Extrapolate the integration point stresses to the nodes using
      // the transposed shape functions.

      ipoint     = ipMpMap_ (ielem,ip);

      material_->update ( stressIp, stiff, strain_(ALL,igpoint) , ipoint );

      ndStress  += matmul ( sfuncs(ALL,ip), stressIp );
      ndWeights += sfuncs(ALL,ip);
    }

    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table.addBlock ( inodes, jcols, ndStress );
  }
}


//-----------------------------------------------------------------------
//   getStrain_
//-----------------------------------------------------------------------


void SaturatedPorousModel::getStrain_

  ( XTable&        table,
    const Vector&  weights )

{
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount = ielems.size         ();
  const int  nodeCount  = shape_->nodeCount   ();
  const int  ipCount    = shape_->ipointCount ();
  const int  strCount   = STRAIN_COUNTS[rank_];

  Matrix     sfuncs     = shape_->getShapeFunctions ();

  Matrix     ndStrain   ( nodeCount, strCount );
  Vector     ndWeights  ( nodeCount );

  IdxVector  inodes     ( nodeCount );
  IdxVector  jcols      ( strCount  );

  int        ipoint;


  // Add the columns for the normal strain components to the table.

  switch ( strCount )
  {
  case 1:

    jcols[0] = table.addColumn ( "e_xx" );

    break;

  case 3:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_xy" );

    break;

  case 4:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_zz" );
    jcols[3] = table.addColumn ( "e_xy" );

    break;  

  case 6:

    jcols[0] = table.addColumn ( "e_xx" );
    jcols[1] = table.addColumn ( "e_yy" );
    jcols[2] = table.addColumn ( "e_zz" );
    jcols[3] = table.addColumn ( "e_xy" );
    jcols[4] = table.addColumn ( "e_yz" );
    jcols[5] = table.addColumn ( "e_xz" );

    break;

  default:

    throw Error (
      JEM_FUNC,
      "unexpected number of strain components: " +
      String ( strCount )
    );
  }

  // Iterate over all elements assigned to this model.

  ipoint = 0;

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    // Get the global element index.

    int  ielem = ielems[ie];

    // Get the element nodes.

    elems_.getElemNodes ( inodes, ielem );

    // Iterate over the integration points.

    ndStrain  = 0.0;
    ndWeights = 0.0;

    for ( int ip = 0; ip < ipCount; ip++, ipoint++ )
    {
      // Extrapolate the integration point strains to the nodes using
      // the transposed shape functions.

      ndStrain  += matmul ( sfuncs(ALL,ip), strain_(ALL,ipoint) );
      ndWeights += sfuncs(ALL,ip);
    }

    // Increment the table weights. When the complete table has been
    // filled, Jive will divide each row in the table by the
    // corresponding table weight. In this way the strain components
    // are automatically averaged over all elements that are attached
    // to a node. The weight vector is initially zero.

    select ( weights, inodes ) += ndWeights;

    // Add the strains to the table.

    table.addBlock ( inodes, jcols, ndStrain );
  }
}

//-----------------------------------------------------------------------
//   getHistory_
//-----------------------------------------------------------------------


void SaturatedPorousModel::getHistory_

  ( XTable&          table,
    const Vector&    weights )
{

}

//-----------------------------------------------------------------------
//   getOutputData_
//-----------------------------------------------------------------------

/** all data at Gauss points are computed, however, only those requested
 *  in the vtk block in the .pro file will be written to vtu files.
 *  Example: vtk.data = "stress_yy | stress_xx" for xx and yy stress 
 *  components.
 */

void SaturatedPorousModel::getOutputData_

  ( Ref<XTable>    table,
    const Vector&  weights,
    const String&  contents,
    const Vector&  state )
          
{
  MChain2    mc2;

  // table filler related stuff

  TbFiller   tbFiller   ( rank_ );

  Slice      iistrain   = tbFiller.announce ( "strain.tensor" ); // iistrain  = [0 1 2 3]
  Slice      iistress   = tbFiller.announce ( "stress.tensor" ); // iistrain  = [0 1 2 3]

  Vector     ipValues   ( tbFiller.typeCount() ); // typeCount() = # of types = 8 in 2D

  Vector     strain     ( ipValues[iistrain]   );
  Vector     stress     ( ipValues[iistress]   );

  // Let TbFiller find out which columns of ndValues to write to 
  // which columns of the table (based on filter in input file)

  IdxVector  i2table;
  IdxVector  jcols;

  tbFiller . setFilter   ( contents ); 
  tbFiller . prepareTable( i2table, jcols, table ); 
  
  IdxVector  ielems     = egroup_.getIndices  ();

  const int  ielemCount  = ielems.size         ();
  const int  nodeCount   = shape_->nodeCount   ();
  const int  ipCount     = shape_->ipointCount ();
  const int  strCount    = STRAIN_COUNTS[rank_];
  const int  dispCount   = nodeCount * rank_;
   
  Matrix     ndValuesOut  ( nodeCount, i2table.size() );
  Matrix     ndValuesOut1 ( nodeCount, i2table.size() );
  Vector     ipValuesOut  ( i2table.size() );

  Matrix      stiff      ( strCount, strCount  );
  Matrix      bd         ( strCount, dispCount );
  Vector      ndWeights  ( nodeCount           ); 
  IdxVector   inodes     ( nodeCount           );


  Properties  params;

  idx_t       ipoint, igpoint = 0;

  Matrix      sfuncs     = shape_->getShapeFunctions ();

  // Iterate over all elements assigned to this model.

  for ( int ie = 0; ie < ielemCount; ie++ )
  {
    int  ielem = ielems[ie];

    elems_.getElemNodes ( inodes, ielem );

    ndValuesOut = 0.;
    ndWeights   = 0.;

    // Loop on integration points 
    
    for ( idx_t ip = 0; ip < ipCount; ip++, igpoint++ )
    {
      ipoint      = ipMpMap_( ielem, ip );
      
 
      strain      = strain_(slice(BEGIN,strCount),igpoint);
      material_->update ( stress, stiff, strain , ipoint );
 
      //funcs = abs ( funcs );

       // apply the filter now, only cols specified by i2table are 
       // written to the table 
      ipValuesOut  = ipValues[i2table]; 

      matmul (ndValuesOut1, sfuncs(ALL,ip), ipValuesOut );


      ndValuesOut += ndValuesOut1;
      ndWeights   += sfuncs(ALL,ip);

    }  // end of loop on integration points
    
    select ( weights, inodes ) += ndWeights;

    // Add the stresses to the table.

    table->addBlock ( inodes, jcols, ndValuesOut );
  }
}

//-----------------------------------------------------------------------
//   checkCommit_
//-----------------------------------------------------------------------

void SaturatedPorousModel::checkCommit_

  ( const Properties&  params )

{
  // System::info() << myName_ << " : check commit ... do nothing!\n";
}

//-----------------------------------------------------------------------
//   setStepSize_
//-----------------------------------------------------------------------

void SaturatedPorousModel::setStepSize_

  ( const Properties&  params )

{
  double       dt;
  params.get ( dt,  XProps::STEP_SIZE   );

  System::debug() << "Setting step size to " << dt << " ("
    << dt / dtime_ * 100 << "\% of previous step size)" << "\n";

  dtime_ = dt;
  
}

//=======================================================================
//   related functions
//=======================================================================

//-----------------------------------------------------------------------
//   newSaturatedPorousModel
//-----------------------------------------------------------------------


static Ref<Model>     newSaturatedPorousModel

  ( const String&       name,
    const Properties&   conf,
    const Properties&   props,
    const Properties&   globdat )

{
  return newInstance<SaturatedPorousModel> ( name, conf, props, globdat );
}


//-----------------------------------------------------------------------
//   declareSaturatedPorousModel
//-----------------------------------------------------------------------


void declareSaturatedPorousModel ()
{
  using jive::model::ModelFactory;

  ModelFactory::declare ( "SaturatedPorous", & newSaturatedPorousModel );
}
