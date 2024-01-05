
#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "Arrays.h"

//========================================================================
//   Basic functions
//========================================================================

double macaulayP  ( double x );
double macaulayN  ( double x );
double heavisideP ( double x );
double heavisideN ( double x );
int    sign       ( double x );

//========================================================================
//   Solve algebraic equations
//========================================================================

// Solving linear equation of one variable
// a*x + b = 0

void                  solveLinearEqn

  ( double&      ans,
    const double a,
    const double b );

// Solving quadratic equation
// a*x^2 + b*x + c = 0

void                  solveQuadEqn

  ( Vector&      ans,
    const double a,
    const double b,
    const double c );

// Solving cubic equation
// a*x^3 + b*x^2 + c*x + d = 0
// only for functions with real roots
// returns the number of roots

idx_t                 solveCubicEqn

  ( Vec3&            ans,
    const double     a,
    const double     b,
    const double     c,
    const double     d);

#endif

