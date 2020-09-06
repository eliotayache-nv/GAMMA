/*
* @Author: eliotayache
* @Date:   2020-06-10 11:18:13
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-06 15:23:32
*/

#include "../fluid.h"
#include "../environment.h"
#include "../err.h"
#include "../constants.h"
#include "../interface.h"
#include "../cell.h"
#include <iostream>
#include <gsl/gsl_roots.h>   // root finding algorithms
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

FluidState::FluidState(){

  for (int i = 0; i < NUM_Q; ++i) prim[i]=0.;
  prim[RHO]=1.;
  prim[PPP]=1.;

  for (int i = 0; i < NUM_Q; ++i) cons[i]=0.;
  prim[DEN]=1.;

  for (int n = 0; n < NUM_D; ++n){
    for (int i = 0; i < NUM_Q; ++i) flux[n][i]=0.;
  }

}

FluidState::~FluidState(){}


void FluidState::prim2cons(){

  double rho = prim[RHO];
  double p   = prim[PPP];
  double u   = 0;
  double uu[NUM_D];
  for (int i = 0; i < NUM_D; ++i){
    uu[i] = prim[UU1+i];
    u += uu[i]*uu[i];
  }
  u = sqrt(u);
  double lfac = sqrt(1+u*u);
  double h   = 1.+p*GAMMA_/(GAMMA_-1.)/rho; // ideal gas EOS (TBC)

  double D   = lfac*rho;
  double tau = D*h*lfac-p-D;
  double ss[NUM_D];
  for (int i = 0; i < NUM_D; ++i) ss[i] = D*h*uu[i];

  cons[DEN] = D;
  cons[TAU] = tau;
  for (int i = 0; i < NUM_D; ++i) cons[SS1+i] = ss[i];

}


void FluidState::state2flux(double r){

  UNUSED(r);

  double p = prim[PPP];
  double u = 0;
  double uu[NUM_D];
  for (int i = 0; i < NUM_D; ++i)
  {
    uu[i] = prim[UU1+i];
    u += uu[i]*uu[i];
  }
  u = sqrt(u);
  double lfac = sqrt(1+u*u);
  double D   = cons[DEN];
  double ss[NUM_D];
  for (int i = 0; i < NUM_D; ++i) ss[i] = cons[SS1+i];

  for (int n = 0; n < NUM_D; ++n){
    flux[n][DEN] = D*uu[n]/lfac;
    for (int i = 0; i < NUM_D; ++i){
      if (i==n) flux[n][SS1+i] = ss[i]*uu[n]/lfac + p;
      else flux[n][SS1+i] = ss[i]*uu[n]/lfac;
    }
    flux[n][TAU] = ss[n] - D*uu[n]/lfac;
    for (int t = 0; t < NUM_T; ++t) flux[n][NUM_C+t] = cons[NUM_C+t]*uu[n]/lfac;
  }

}

/**
 *
 * PRIMITIVE VARIABLES RECONSTRUCTION:
 *
 */

struct f_params{
  double D,S,E,gamma;
};

// This is the function we need to set to zero:
static double f(double p, void *params){

  double lfac; 
  // parsing parameters
  struct f_params *par = (struct f_params *) params;
  double  D = par->D;
  double  S = par->S;
  double  E = par->E;
  double  gamma = par->gamma;

  lfac = 1. / sqrt(1. - (S * S) / ((E + p) * (E + p)));

  return (E + p - D * lfac - (gamma * p * lfac * lfac) / (gamma - 1.));
    // Mignone (2006) eq. 5

}

void FluidState::cons2prim(double pin){

  double  lfac,E,S,D;  // (computed from conserved variables)
  int     status;
  int     iter = 0, max_iter = 1000000;
  double  r;
  double  p_lo, p_hi,p;
  struct  f_params                params;
  const   gsl_root_fsolver_type   *T;
  gsl_root_fsolver                *s;
  gsl_function                    F;

  F.function = &f;
  F.params = &params;

  // setting initial parameters
  D = cons[DEN];
  S = 0.;
  for (int i = 0; i < NUM_D; ++i) S += cons[SS1+i]*cons[SS1+i];
  S = sqrt(S);
  E = cons[TAU]+cons[DEN];

  params.D = D;
  params.S = S;
  params.E = E;
  params.gamma = GAMMA_;

  if (pin != 0) p = pin;
  else p = prim[PPP];

  // Looking for pressure only if current one doesn't work;
  double f_init = f(p, &params);

  if (fabs(f_init) > 1.e-14) {

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);

    // setting boundaries
    // f(p_lo) has to be positive. No solution otherwise
    // hence, f(p_hi) has to be negative

    p_lo = fmax(0., (1. + 1e-13) * fabs(S) - E);

    if (f_init < 0){
      p_hi = 1.*p;
    }
    else {
      int i = 0;
      p_hi = p;
      while (f(p_hi, &params) > 0) {
        i++;
        p_hi *= 10;
      }
    }
    gsl_root_fsolver_set (s, &F, p_lo, p_hi);

    // Iterating for root finding:            
    do {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      p_lo = gsl_root_fsolver_x_lower (s);
      p_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (p_lo, p_hi,
                       0, 1.e-13);
    }                
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);
    p = r;
  }

  lfac = 1./sqrt(1.-(S*S)/((E+p)*(E+p)));

  prim[PPP] = p;
  prim[RHO] = D/lfac;

  // printf("%le\n", (S*S)/((E+p)*(E+p)));

  if (S == 0.)    
    for (int i = 0; i < NUM_D; ++i) prim[UU1+i] = 0;
  else
  {
    double v = sqrt(1.- 1./(lfac*lfac));
    for (int i = 0; i < NUM_D; ++i) prim[UU1+i] = lfac*v*(cons[SS1+i]/S);
    // Mignone (2006) eq. 3
  }

  prim2cons(); // for consistency
  
}

void Interface::wavespeedEstimates(){

  int u = UU1+dim;
  // if      (dim == x_) { u = UU1; }
  // else if (dim == y_) { u = UU2; }
  // else                { u = UU3; }

  double lfacL = SL.lfac();
  double lfacR = SR.lfac();
  double vL = SL.prim[u]/lfacL;
  double vR = SR.prim[u]/lfacR;

  // Computing speed of sound parameter sigmaS:
  double cSL = sqrt(GAMMA_ * SL.prim[PPP] 
    / (SL.prim[RHO] + SL.prim[PPP]*GAMMA_/(GAMMA_-1.)));
  double cSR = sqrt(GAMMA_ * SR.prim[PPP] 
    / (SR.prim[RHO] + SR.prim[PPP]*GAMMA_/(GAMMA_-1.)));
    // Mignone (2006) eq. 4
  double sgSL = cSL * cSL / (lfacL * lfacL * (1. - cSL * cSL));
  double sgSR = cSR * cSR / (lfacR * lfacR * (1. - cSR * cSR));
    // Mignone (2006) eq. 22-23s

  // Retreiving L and R speeds:
  double  l1, l2;       // 2 intermediates
  l1 = (vR - sqrt(sgSR * (1. - vR * vR + sgSR))) / (1. + sgSR);
  l2 = (vL - sqrt(sgSL * (1. - vL * vL + sgSL))) / (1. + sgSL);
  lL = fmin(l1,l2);

  l1 = (vR + sqrt(sgSR * (1. - vR * vR + sgSR))) / (1. + sgSR);
  l2 = (vL + sqrt(sgSL * (1. - vL * vL + sgSL))) / (1. + sgSL);
  lR = fmax(l1,l2);

}

void Interface::computeLambda(){

  int u = UU1+dim;
  int s = SS1+dim;
  // if      (dim == x_) { u = UU1; s = SS1; }
  // else if (dim == y_) { u = UU2; s = SS2; }
  // else                { u = UU3; s = SS3; }

  wavespeedEstimates();
  
  double lfacL = SL.lfac();
  double lfacR = SR.lfac();
  double vL = SL.prim[u]/lfacL;
  double vR = SR.prim[u]/lfacR;
  double mL = SL.cons[s];
  double mR = SR.cons[s];
  double pL = SL.prim[PPP];
  double pR = SR.prim[PPP];
  double EL = SL.cons[TAU]+SL.cons[DEN];
  double ER = SR.cons[TAU]+SR.cons[DEN];

  // Setting up temporary variables:
  double AL = lL * EL - mL;
  double AR = lR * ER - mR;
  double BL = mL * (lL - vL) - pL;
  double BR = mR * (lR - vR) - pR;
    // Mignone (2006) eq. 17

  // Coefficients for lambdaStar equation:
  double FhllE = (lL * AR - lR * AL) / (lR - lL);
  double Ehll  = (AR - AL) / (lR - lL);
  double Fhllm = (lL * BR - lR * BL) / (lR - lL);
  double mhll  = (BR - BL) / (lR - lL);

  if (FhllE == 0){
    lS = mhll / (Ehll + Fhllm);
    return;
  }

  double delta = pow(Ehll + Fhllm,2.) - 4. * FhllE * mhll;
  lS = ((Ehll + Fhllm) - sqrt(delta)) / (2. * FhllE);

}


FluidState Interface::starState(FluidState Sin, double lbda){

  int u,s, uu[2], ss[2];
  if (dim == x_) {
    u = UU1 + x_;
    s = SS1 + x_;
    uu[0] = UU1 + y_; uu[1] = UU1 + z_;
    ss[0] = SS1 + y_; ss[1] = SS1 + z_;
  }
  else if (dim == y_) {
    u = UU1 + y_;
    s = SS1 + y_;
    uu[0] = UU1 + x_; uu[1] = UU1 + z_;
    ss[0] = SS1 + x_; ss[1] = SS1 + z_;
  }
  else {
    u = UU1 + z_;
    s = SS1 + z_;
    uu[0] = UU1 + x_; uu[1] = UU1 + y_;
    ss[0] = SS1 + x_; ss[1] = SS1 + y_;
  }

  double lfac = Sin.lfac();
  double v = Sin.prim[u]/lfac;
  double m = Sin.cons[s];
  double p = Sin.prim[PPP];
  double D = Sin.cons[DEN];
  double E = Sin.cons[TAU]+Sin.cons[DEN];

  // Setting up temporary variables:
  double A = lbda * E - m;
  double B = m * (lbda-v) - p;
    // Mignone (2006) eq. 17
  FluidState Sout;
  double pS = (A*lS - B) / (1. - lbda*lS);
    // Mignone (2006) eq. 17
  double k = (lbda-v) / (lbda - lS);
  Sout.cons[DEN] = D * k;
  Sout.cons[s]   = (m * (lbda-v) + pS - p) / (lbda-lS); 
  for (int n = 0; n < NUM_D-1; ++n){ Sout.cons[ss[n]] = Sin.cons[ss[n]] * k; }
  Sout.cons[TAU] = (E * (lbda-v) + pS * lS - p * v) / (lbda - lS) - D*k;
      // Mignone (2006) eq. 16 (-D)

  return Sout;
  // no need to compute the flux, will be done in the solver function

}


void Cell::sourceTerms(double dt){

  return;

}

