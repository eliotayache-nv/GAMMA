/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-11-28 11:37:29
*/

#include "../environment.h"
#include "../grid.h"
#include <gsl/gsl_roots.h>   // root finding algorithms
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>


void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 100;
  par->ncell[y_] = 100;
  par->nmax      = 1000;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){
  // Careful! When switching moving coordinate, you need to set MV as a function of i
  // and F1 as a functino of j

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) 2.*(i+0.5)/ncell[x_] - 1.;
      c->G.dx[x_] =          2./ncell[x_];
      c->G.x[y_]  = (double) 2.*(j+0.5)/ncell[y_] - 1.;
      c->G.dx[y_] =          2./ncell[y_];

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (int j = 0; j < ncell[x_]; ++j){
    for (int i = 0; i < ncell[y_]; ++i){
      Cell *c = &Cinit[j][i];

      double x = c->G.x[x_];
      double y = c->G.x[y_];

      double rho1 = 1.;
      double rho2 = 0.1;
      double p0 = 1;

      c->S.prim[RHO] = rho1;
      c->S.prim[VV1] = 0.;
      c->S.prim[VV2] = 0.;
      c->S.prim[PPP] = p0;
      c->S.prim[TR1] = 2.;

      if (fabs(y) < 0.2 and x<-0.5){
        c->S.prim[RHO] = rho1;
        c->S.prim[VV1] = 0.99999;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = p0;
        c->S.prim[TR1] = 1.;

      }

      // c->S.prim[VV1] = 0.02 * (double) gsl_rng_get (r) / (double) gsl_rng_max(r) - 0.01;

    }
  }

  gsl_rng_free (r);
  return 0;

}


void Grid::userKinematics(int it, double t){

  // setting lower and higher i boundary interface velocities to zero
  double vIn     = 0.3;
  double vOut    = 1.05;     
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n < ngst; ++n){
      int    iL = n;
      int    iR = ntrack[j]-2-n;
      Itot[j][iL].v = vIn;
      Itot[j][iR].v = vOut;
    }
  }

}


void Cell::userSourceTerms(double dt){

  // // sources have to be expressed in terms of conserved variables
  // double y   = G.x[y_]; 
  // double dV  = G.dV;
  // double rho = S.prim[RHO];
  // double fg = -rho*g*dV*dt;

  // S.cons[SS1] += fg;

}


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell c = Ctot[j][i];

  double split_dl = 0.05;
  double merge_dl = 0.005;
    // careful, dx != dl
  double x = c.G.x[x_];
  double dx = c.G.dx[x_];
  double x1 = Ctot[j-1][iRbnd[j-1]-1].G.x[x_];
  double dx1 = Ctot[j-1][iRbnd[j-1]-1].G.dx[x_];

  if (c.G.dx[MV] > split_dl) {
    printf("split %d %d %d %d %le %le \n", j, i, iLbnd[j], iRbnd[j], x+dx/2., x1+dx1/2.);
    return(split_);
  }
  if (c.G.dx[MV] < merge_dl) {
    return(merge_);
  }
  return(skip_);

}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value
  
  UNUSED(*res);

}


void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  UNUSED(*rho);
  UNUSED(*p);

  return;

}














