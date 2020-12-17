/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-12-16 19:18:05
*/

#include "../../environment.h"
#include "../../grid.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 200;
  par->ncell[y_] = 100;
  par->nmax      = 400;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){
  // Careful! When switching moving coordinate, you need to set MV as a function of i
  // and F1 as a functino of j

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) 2*(i+0.5)/ncell[x_];
      c->G.dx[x_] =          2./ncell[x_];
      c->G.x[y_]  = (double) (j+0.5)/ncell[y_] - 0.5;
      c->G.dx[y_] =          1./ncell[y_];

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double x = c->G.x[x_];
      double y = c->G.x[y_];
      if (fabs(y) < 0.2 and x < 0.3){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.9;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 1.;
        c->S.cons[NUM_C] = 1.;
      }
      else {
        c->S.prim[RHO] = 1.;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 1.;
        c->S.cons[NUM_C] = 2.;
      }
    }

  }

  return 0;

}


void Grid::userKinematics(int it, double t){

  UNUSED(it);
  UNUSED(t);

  // setting lower and higher i boundary interface velocities to zero
  double vIn     = 0.;     // can't be lower than 1 for algo to work
  double vOut    = 0.;     
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

  UNUSED(dt);

}


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell c = Ctot[j][i];

  double split_dl = 0.02;
  double merge_dl = 0.005;
    // careful, dx != dl

  if (c.G.dx[MV] > split_dl) {
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














