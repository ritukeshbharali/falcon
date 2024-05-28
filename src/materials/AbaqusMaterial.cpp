#include <jem/base/System.h>
#include <jem/base/limits.h>
#include <jem/base/Float.h>
#include <jem/base/Error.h>
#include <jem/base/PrecheckException.h>
#include <jem/base/Array.h>
#include <jem/util/Properties.h>
#include <jem/numeric/algebra/utilities.h>
#include <jem/numeric/utilities.h>
#include <jem/io/NumberFormat.h>
#include <algorithm>
#include <dlfcn.h>
#include <cstdlib>
#include <iostream>

#include "util/Constants.h"
#include "AbaqusMaterial.h"

using namespace jem;
using namespace jem::io;

using jem::util::Properties;
using jive::Vector;


//=======================================================================
//   class AbaqusMaterial
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char*  AbaqusMaterial::UMAT_MATPROPS     = "matProps";
const char*  AbaqusMaterial::UMAT_NSTATES      = "nStates";
const char*  AbaqusMaterial::STORE_OLD_STRAIN  = "storeOldStrain";

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


AbaqusMaterial::AbaqusMaterial

  ( int rank, const Properties& globdat )
    : Material ( rank, globdat ), nStates_ ( 0 )
{
  // Attempt to load the shared library 'libumat.so'
  // RLTD_NOW resolves all symbols upon loading.
  void* handle = dlopen(UMAT_SHARED_LIB, RTLD_NOW);

  // Check if 'handle' is a valid pointer
  if ( handle == nullptr )
  {
    std::cerr << "Error loading shared library: "
              << dlerror()
              << std::endl;
  }

  // Function pointer to UMAT
  umat = reinterpret_cast<umatFunc>
              ( dlsym( handle, "umat_" ) );

  // Check function pointer validity
  if ( umat == nullptr )
  {
    std::cerr << "Error locating UMAT function: "
              << dlerror()
              << std::endl;
  }

  // Abaqus strain vector is ordered differently from
  // that of falcon. Therefore, perm_ is required to
  // store the array permutations.
  // Abaqus: [xx, yy, zz, xy, xz, yz]
  // falcon: [xx, yy, zz, xy, yz, xz]
  perm_.resize(6);
  perm_ = {0,1,2,3,5,4};

  // Set history allocation flag to false
  allocd_ = false;
}


AbaqusMaterial::~AbaqusMaterial ()
{} 


//-----------------------------------------------------------------------
//   configure
//-----------------------------------------------------------------------


void AbaqusMaterial::configure ( const Properties& props,
                                 const Properties& globdat )
{
  props.get ( matProps_, UMAT_MATPROPS    );
  props.get ( nStates_,  UMAT_NSTATES     );
  props.get ( storeOld_, STORE_OLD_STRAIN );
}


//-----------------------------------------------------------------------
//   getConfig
//-----------------------------------------------------------------------


void AbaqusMaterial::getConfig ( const Properties& conf,
                                 const Properties& globdat ) const
{
  conf.set ( UMAT_MATPROPS,    matProps_ );
  conf.set ( UMAT_NSTATES,     nStates_  );
  conf.set ( STORE_OLD_STRAIN, storeOld_ );
}


//-----------------------------------------------------------------------
//   updateConfig
//-----------------------------------------------------------------------

void AbaqusMaterial::updateConfig ()

{
  // do nothing!
}

//-----------------------------------------------------------------------
//   update
//-----------------------------------------------------------------------

void  AbaqusMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      idx_t                 ipoint )
{
  Vector dstrain; dstrain.resize(strain.size()); dstrain = 0.0;

  if ( allocd_ )
  {
    dstrain = strain - preStrain_[ipoint].v;
  }

  update ( stress, stiff, strain, dstrain, ipoint );
}

// ----------------------------------------------------------------
//  update (overloaded version)
// ----------------------------------------------------------------

void  AbaqusMaterial::update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const Vector&         dstrain,
      int                   ipoint )
{
  // Some constants

  const int strCount   = STRAIN_COUNTS[rank_];
  const int matCount   = matProps_.size();

  // Initialize UMAT variables

  double STRESS[strCount] = {0.0};            // Zero vector
  double STATEV[nStates_] = {0.0};            // Zero vector
  double DDSDDE[strCount][strCount] = {0.0};  // Zero matrix
  double SSE;                                 // Not used
  double SPD;                                 // Not used
  double SCD;                                 // Not used
  double RPL;                                 // Not used
  double DDSDDT;                              // Not used
  double DRPLDE;                              // Not used
  double DRPLDT;                              // Not used
  double STRAN[strCount]  = {0.0};            // Zero vector
  double DSTRAN[strCount] = {0.0};            // Zero vector
  double TIME;                                // Not used
  double DTIME;                               // Not used
  double TEMP;                                // Not used
  double DTEMP;                               // Not used
  double PREDEF;                              // Not used
  double DPRED;                               // Not used
  char   CMNAME;                              // Not used 
  int    NDI   = 3;                           // 3 Direct stress
  int    NSHR  = strCount - NDI;              // Shear stress
  int    NTENS = strCount;                    // (= NDI + NSHR)        
  int    NSTATV;                              // Not used int        
  double PROPS[matCount] = {0.0};             // Initialize mat props to zero
  int    NPROPS = matCount;                   // # material properties
  double COORDS;                              // Not used       
  double DROT;                                // Not used       
  double PNEWDT;                              // Not used       
  double CELENT;                              // Not used       
  double DFGRD0;                              // Not used       
  double DFGRD1;                              // Not used       
  int    NOEL;                                // Not used       
  int    NPT;                                 // Not used       
  int    LAYER;                               // Not used       
  int    KSPT;                                // Not used       
  int    KSTEP;                               // Not used       
  int    KIN;                                 // Not used

  // Fill material properties
  for ( int i = 0; i < matCount; i++ )
  {
    PROPS[i] = matProps_[i];
  }

  // Fill strain and strain increment
  for ( int i = 0; i < strCount; i++ )
  {
    STRAN[i]  = strain [perm_[i]];
    DSTRAN[i] = dstrain[perm_[i]];
  }

  // Fill state variables
  if ( allocd_ )
  {
    for ( int i = 0; i < nStates_; i++ )
    {
      STATEV[i] = preState_[ipoint].v[i];
    }
  }

  // Call Abaqus UMAT function
  umat(&STRESS[0], &STATEV[0], &DDSDDE[0][0], &SSE, &SPD,
       &SCD, &RPL, &DDSDDT, &DRPLDE, &DRPLDT, &STRAN[0], 
       &DSTRAN[0], &TIME, &DTIME, &TEMP, &DTEMP, &PREDEF, 
       &DPRED, &CMNAME, &NDI, &NSHR, &NTENS, &NSTATV, 
       &PROPS[0], &NPROPS, &COORDS, &DROT, &PNEWDT, &CELENT,
       &DFGRD0, &DFGRD1, &NOEL, &NPT, &LAYER, &KSPT, &KSTEP,
       &KIN);

  // Transfer Abaqus output to Jive variables (stress, stiff)
  for ( int i = 0; i < strCount; i++ )
  {
    stress[i] = STRESS[perm_[i]];

    for ( int j = 0; j < strCount; j++ )
    {
      stiff(i,j) = DDSDDE[perm_[i]][perm_[j]];
    }
  }

  // Update state variables
  if ( allocd_ )
  {
    for ( int i = 0; i < nStates_; i++ )
    {
      newState_[ipoint].v[i] = STATEV[i];
    }
  }
}


//-----------------------------------------------------------------------
//   allocPoints
//-----------------------------------------------------------------------

void AbaqusMaterial::allocPoints ( int count )
{
  preState_.reserve ( count );
  newState_.reserve ( count );

  for ( idx_t ip = 0; ip < count; ++ip )
  {
    preState_.pushBack ( Hist_( nStates_ ) );
    newState_.pushBack ( Hist_( nStates_ ) );
  }

  latestState_ = &preState_;

  if ( storeOld_ )
  {
    preStrain_.reserve ( count );
    newStrain_.reserve ( count );

    for ( idx_t ip = 0; ip < count; ++ip )
    {
      preStrain_.pushBack ( Hist_( STRAIN_COUNTS[rank_] ) );
      newStrain_.pushBack ( Hist_( STRAIN_COUNTS[rank_] ) );
    }

    latestStrain_ = &preStrain_;
  }

  // Set history allocation flag to true
  allocd_ = true;
}


// --------------------------------------------------------------------
//  commit
// --------------------------------------------------------------------

void  AbaqusMaterial::commit()
{
  newState_.swap  ( preState_ );
  latestState_ = &( preState_ );

  if ( storeOld_ )
  {
    newStrain_.swap  ( preStrain_ );
    latestStrain_ = &( preStrain_ );
  }
}


//-----------------------------------------------------------------------
//  clone
//-----------------------------------------------------------------------

Ref<Material> AbaqusMaterial::clone ( ) const

{
  return newInstance<AbaqusMaterial> ( *this );
}


//-----------------------------------------------------------------------
//  Hist_ constructor
//-----------------------------------------------------------------------

AbaqusMaterial::Hist_::Hist_ ( int n )
{
  v.resize ( n ); v = 0.0;
}