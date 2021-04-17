/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-04-11 21:04:32
*/

#include "../../environment.h"
#include "../../grid.h"
#include "../../simu.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 300;
  par->ncell[y_] = 300;
  par->nmax      = 600;    // max number of cells in MV direction
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
      c->G.x[y_]  = (double) 2.*(j+0.5)/ncell[y_] - 1;
      // c->G.x[y_]  = (double) 2.*(j+0.5)/ncell[y_] - 1.000001;
      c->G.dx[y_] =          2./ncell[y_];

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
      if (x>0 and y>0){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 0.01;
      }
      if (x<0 and 0<=y){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.99;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 1.;
      }
      if (x<=0 and y<=0){
        c->S.prim[RHO] = 0.5;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 1.;
      }
      if (x>0 and 0>=y){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0;
        c->S.prim[VV2] = 0.99;
        c->S.prim[PPP] = 1.;
      }
      // if (y > x and y >= -x){
      //   c->S.prim[RHO] = 0.1;
      //   c->S.prim[VV1] = 0.;
      //   c->S.prim[VV2] = 0.;
      //   c->S.prim[PPP] = 0.01;
      // }
      // if (y < -x and y > x){
      //   c->S.prim[RHO] = 0.1;
      //   c->S.prim[VV1] = 0.7;
      //   c->S.prim[VV2] = 0.7;
      //   c->S.prim[PPP] = 1.;
      // }
      // if (y <= x and y < -x){
      //   c->S.prim[RHO] = 0.5;
      //   c->S.prim[VV1] = 0.;
      //   c->S.prim[VV2] = 0.;
      //   c->S.prim[PPP] = 1.;
      // }
      // if (y >= -x and y <= x){
      //   c->S.prim[RHO] = 0.1;
      //   c->S.prim[VV1] = -0.7;
      //   c->S.prim[VV2] = 0.7;
      //   c->S.prim[PPP] = 1.;
      // }
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

  double split_ar = 3;
  double merge_ar = 0.2;
  double dx = c.G.dx[x_];
  double dy = c.G.dx[y_];
  double ar = dx/dy;

  if (ar > split_ar) {
    return(split_);
  }
  if (ar < merge_ar) {
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




void Simu::dataDump(){

  // if (it%1 == 0){ grid.printCols(it, t); }
  if (it%100 == 0){ grid.printCols(it, t); }

}

void Simu::runInfo(){

  // if ((worldrank == 0) and (it%1 == 0)){ printf("it: %ld time: %le\n", it, t);}
  if ((worldrank == 0) and (it%100 == 0)){ printf("it: %ld time: %le\n", it, t);}

}

void Simu::evalEnd(){

  // if (it > 300){ stop = true; }
  if (t > 0.8){ grid.printCols(it, t); stop = true; } // 3.33e8 BOXFIT simu

}















