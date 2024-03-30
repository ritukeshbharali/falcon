/*
 * 
 *  Copyright (C) 2007 TU Delft. All rights reserved.
 *  
 *  This class implemens the isotropic elastic-damage material
 *  This represents the material at a point in space.
 *  It is implemented in such a way that can be used for any
 *  discretisation methods, say finite elements, EFG and so on.
 *  
 *  Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *
 */

/** @file DamageMaterial.h
 *  @brief Isotropic damage material model.
 * 
 *  Original Author: V.P. Nguyen, V.P.Nguyen@tudelft.nl
 *  Date: 30 November 2007
 *  
 *  Updates (when, what and who)
 *     - [21 February 2024] Integrated old damage laws,
 *     linear1, linear2, and new polynomial and multilinear.
 *     Removed Rankine type equivalent strain measure. (RB)
 */

#ifndef DAMAGEMATERIAL_H
#define DAMAGEMATERIAL_H

#include <jive/Array.h>
#include <jem/base/String.h>
#include <jem/util/Flex.h>
#include <vector>

#include "HookeMaterial.h"
#include "VonMisesStrain.h"


using jem::String;
using jem::util::Flex;
using jive::Vector;
using jive::IntVector;
using std::vector;


enum SofteningLaw {
  Perfect,
  Linear1,
  Linear2,
  Exponential3Params,
  Exponential2Params,
  Exponential2ParamsReg,
  Hyperbolic,
  Power,
  ExpoEnergy,
  Polynomial,
  Multilinear
};

enum EquiStrainDef {
  Mazars,
  VonMises,
  StrainEnergy
};


// =======================================================================
//  class DamageMaterial
// =======================================================================

/** @brief 
 *  The DamageMaterial class implements an isotropic damage material.
 *  stress = (1-omega) D:strain
 * 
 *  Original Author: Vinh Phu Nguyen (Aug, 2007, TU Delft)
 *  See ofeFRAC (https://github.com/vinhphunguyen/ofeFRAC)
 * 
 *  Contributions from RB's M.Sc thesis work (Feb 2017, TU Delft) added.
 *  The separate damage material laws were merged into Nguyen's code.
 *      - LinSoftIsoDamage (Linear1)
 *      - LinIsoDamage     (Linear2)
 *      - ExpIsoDamage     (Exponential1)
 * 
 *  New softening laws added:
 *      - Polynomial       (COMSOL Manual)
 *      - Multilinear      (COMSOL Manual)
 *  
 */


class DamageMaterial : public Material
{
 public:

  typedef  DamageMaterial Self;

  static const char*      SOFTENING_PROP;
  static const char*      EQUISTRAIN_PROP;

  static const char*      KAPPAI_PROP;
  static const char*      KAPPAP_PROP;
  static const char*      KAPPAC_PROP;
  static const char*      ALPHA_PROP;
  static const char*      BETA_PROP;
  static const char*      ETA_PROP;
  static const char*      B_PROP;
  static const char*      LENGTH_PROP;

  static const char*      TENSILE_PROP;
  static const char*      FRACTURE_ENERGY_PROP;

  static const char*      CRACK_WIDTH_PROP;
  static const char*      REMOVE_DAMAGE_PROP;

  static const char*      MAZARS_EQUI_STRAIN;
  static const char*      MISES_EQUI_STRAIN;
  static const char*      RANKINE_EQUI_STRAIN;
  static const char*      ENERGY_STRAIN;

  static const char*      PERFECT_SOFTENING;
  static const char*      LINEAR1_SOFTENING;
  static const char*      LINEAR2_SOFTENING;
  static const char*      EXPONENT1_SOFTENING; 
  static const char*      EXPONENT2_SOFTENING;
  static const char*      EXPONENT3_SOFTENING;
  static const char*      POWER_SOFTENING;
  static const char*      HYPERBOLIC_SOFTENING;
  static const char*      EXPOENERGY_SOFTENING;
  static const char*      POLYNOMIAL_SOFTENING;
  static const char*      MULTILINEAR_SOFTENING;

  static const double     CRITICAL_DAMAGE;

  static vector<String>   equiStrainDefs;
  static vector<String>   softeningLawDefs;

  explicit                DamageMaterial

    ( int                   rank,
      const Properties&     globdat );

  virtual void            configure

    ( const Properties&     props,
      const Properties&     globdat );

  virtual void            getConfig

    ( const Properties&     conf,
      const Properties&     globdat ) const;

  virtual void            updateConfig (); 

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      const idx_t           ipoint );   

  // the same as above but overloaded for mesh adjusted damage model
  // he is the element characteristic element length 

  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain,
      int                   ipoint,
      double                he );
  
  virtual void            update

    ( Vector&               stress,
      Matrix&               stiff,
      const Vector&         strain0,
      const Vector&         dstrain,
      int                   ipoint,
      double                he );

  virtual void            allocPoints

    ( int                   ipCount );

  virtual Ref<Material>   clone  () const;

  virtual void            commit ();

  // ----------------------------------------------------------
  //  GETTERS
  // ----------------------------------------------------------

  // return internal variable kappa

  inline double           giveHistory      ( int ipoint ) const;
  inline Flex<double>     giveHistory      (            ) const;

  // return the damage variable

  inline double           giveOmega        ( int ipoint ) const;
  inline Vector           giveOmega        ()             const;

  // return true/false based on damage evolution

  int                     isLoading        ( int ipoint ) const;

  // return true/false if material point is fully damaged

  bool                    isFullyDamaged   ( int ipoint ) const;

  // return the length scale

  inline double           giveLengthScale  ()             const;

  // return the elastic material stiffness matrix

  Matrix                  giveElasticMatrix ()            const;

  // compute the equivalent strain given a strain vector

  double                  getEquiStrain

    ( const Vector& strain );

  // compute derivatives of equivalent strain w.r.t strain
  // vector

  Vector                  getdEpsBardEps

    ( const Vector& strain );  

  // compute the derivative of Omega w.r.t kappa

  double                  getdOmegadKappa

    ( double kappa )                                      const;

  double                  getdOmegadKappa

    ( double kappa, double he )                           const;
  
 
 protected:

  virtual                ~DamageMaterial   ();

 private:

  // compute the damage variable omega from kappa

  double                  damageEvolution_

    ( double kappa )                                      const;

  // compute the damage variable omega from kappa
  // and mesh adjusted size parameter 'he'

  double                  damageEvolution_

    ( double kappa,
      double he )                    const;
  

  // ---------------------------------------------------------------
  // compute derivatives of Mazars equivalent strain w.r.t strain
  // ---------------------------------------------------------------

  Vector                  getdEquidEpsMazars_

    ( const Vector&   strain );

  
 private:

  SofteningLaw            softening_;    //! Softening law 
  EquiStrainDef           equiStrain_;   //! Equivalent strain

  double                  kappaI_;       //! Damage initiation threshold
  double                  kappaP_;       //! Multilinear softening param
  double                  kappaC_;       //! Softening modulus
  double                  alpha_;        //! Exponential softening param
  double                  beta_;         //! Exponential softening param
  double                  b_;            //! Hyperbolic softening param
  double                  s_;            //! ExpoEnergy softening param 
  double                  eta_;          //! Ratio fc/ft
  double                  threshold_;    //! isFullyDamaged threshold

  double                  c_;            //! length-scale squared
  double                  ft_;           //! tensile strength
  double                  gf_;           //! fracture energy

  double                  lambda_;       //! crack band width 

  /**
   *  Struct to store internal variables (equivalent strain, loading)
   */ 

  struct                  hist_
  {
    Flex<double>            eqveps ;     //! non-local equivalent strain
    Flex<int>               loading;     //! loading/unloading
  };

  hist_                   preHist_;      //! previous load step
  hist_                   newHist_;      //! current iteration

  Ref<HookeMaterial>      elasticMat_;   //! pointer to HookeMaterial
  Matrix                  elasticMod_;   //! elastic material stiffness

  VonMisesStrain          vonMisesEqv_;  //! instance of VonMisesStrain

};


// -------------------------------------------------------------------
//  giveHistory ( return kappa for ipoint )
// -------------------------------------------------------------------

inline double  DamageMaterial::giveHistory ( int ipoint ) const
{
  return newHist_.eqveps[ipoint];
}


// -------------------------------------------------------------------
//  giveHistory ( returns kappa for all ipoints )
// -------------------------------------------------------------------

inline Flex<double>  DamageMaterial::giveHistory ( ) const
{
  return newHist_.eqveps;
}


// -------------------------------------------------------------------
//  giveOmega ( returns damage for ipoint )
// -------------------------------------------------------------------

inline double  DamageMaterial::giveOmega ( int ipoint )  const
{
  return damageEvolution_ ( newHist_.eqveps[ipoint] ); 
}


// -------------------------------------------------------------------
//  giveOmega ( returns damage for all ipoints )
// -------------------------------------------------------------------

inline Vector  DamageMaterial::giveOmega()  const
{
  const int s =  preHist_.eqveps.size ( );

  Vector omega ( s );

  for ( int i = 0; i < s; i++ )
  {
    omega[i] =  damageEvolution_ ( newHist_.eqveps[i] ); 
  }

  return omega;
}


// -------------------------------------------------------------------
//  isLoading
// -------------------------------------------------------------------

inline int    DamageMaterial::isLoading ( int ipoint ) const
{
  return newHist_.loading[ipoint];
}


// -------------------------------------------------------------------
//  giveLengthScale
// -------------------------------------------------------------------

inline double  DamageMaterial::giveLengthScale   ( ) const
{
  return c_;
}


#endif


// ======================================================================
//   Softening laws
// ======================================================================

// ----------------------------------------------------------------------
//  perfect softening 
// ----------------------------------------------------------------------

double                  perfectSoftening

  ( double Ki,
    double K );

double                  getDerivPerfectSoftening

    ( double Ki,
      double K );

// ----------------------------------------------------------------------
//  linear softening (type I)
// ----------------------------------------------------------------------

double                  linear1Softening

  ( double Ki,
    double Kc,
    double K );

double                  getDerivLinear1Softening

    ( double Ki,
      double Kc,
      double K );

// ----------------------------------------------------------------------
//  linear softening (type II)
// ----------------------------------------------------------------------

double                  linear2Softening

  ( double Ki,
    double Kc,
    double K );

double                  getDerivLinear2Softening

    ( double Ki,
      double Kc,
      double K );

// ----------------------------------------------------------------------
//  exponential softening1 ( Ki, alpha, beta, K) 
// ----------------------------------------------------------------------

double                  exponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    const double K );

double                  getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double beta,
    double K );

// overloaded for mesh adjusted softenting 

double                  exponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K, 
    double he );

double                  getDerivExponentialSoftening1

  ( double Ki,
    double alpha,
    double gf, double ft,
    double K,
    double he );

// ----------------------------------------------------------------------
//  exponential softening2 
// ----------------------------------------------------------------------2

double                  exponentialSoftening2

  ( double Ki,
    double Kc,
    double K );

double                  getDerivExponentialSoftening2

    ( double Ki,
      double Kc,
      double K );

// overloaded for mesh adjusted softening modulus

double                  exponentialSoftening2

  ( double Ki, double Kc,
    double K, double lamda,
    double he );

double                  getDerivExponentialSoftening2

    ( double Ki,double Kc,
      double K, double lamda,
      double he );

// Jirasek regularised exponential softening law

double                  exponentialSoftening3

  ( double Ki, double K,
    double gf, double ft, 
    double he );

double                  getDerivExponentialSoftening3

   ( double Ki, double K,
     double gf, double ft, 
     double he );


// ----------------------------------------------------------------------
//  power softening 
// ----------------------------------------------------------------------

double                  powerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K );

double                  getDerivPowerSoftening

  ( double Ki,
    double Kc,
    double alpha,
    double beta,
    double K );

// ----------------------------------------------------------------------
//  hyperbolic softening 
// ----------------------------------------------------------------------

double                  hyperbolicSoftening

  ( double Ki,
    double b,
    double K );

double                  getDerivHyperbolicSoftening

  ( double Ki,
    double b,
    double K );


// ----------------------------------------------------------------------
//  exponential strain energy softening 
// ----------------------------------------------------------------------

double                     expoEnergySoftening

 ( double kappa, double kappaI,
   double kappaC, double s );

double                     getDerivExpoEnergySoftening

 ( double kappa, double kappaI,
   double kappaC, double s );


// ----------------------------------------------------------------------
//  polynomial softening (COMSOL Manual)
// ----------------------------------------------------------------------

double                  polynomialSoftening

  ( const double Ki,
    const double Kc,
    const double K );

double                  getDerivPolynomialSoftening

    ( const double Ki,
      const double Kc,
      const double K );


// ----------------------------------------------------------------------
//  multilinear softening (COMSOL Manual)
// ----------------------------------------------------------------------

double                  multilinearSoftening

  ( const double Ki,
    const double Kp,
    const double Kc,
    const double K );

double                  getDerivMultilinearSoftening

    ( const double Ki,
      const double Kp,
      const double Kc,
      const double K );
