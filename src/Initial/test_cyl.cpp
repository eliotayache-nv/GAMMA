/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-11 23:00:15
*/

#include "../environment.h"
#include "../grid.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 100;
  par->ncell[y_] = 100;
  par->nmax      = 1000;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){
  // Careful! the moving dimension must evolve with j
  // The fixed dimension must evolve with i

  for (int j = 0; j < ncell[r_]; ++j){
    for (int i = 0; i < ncell[zcyl_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[r_]     = (double) 2.*(j+0.5)/ncell[r_];
      c->G.dx[r_]    =          2./ncell[r_];
      c->G.x[zcyl_]  = (double) 2.*(i+0.5)/ncell[zcyl_] - 1.;
      c->G.dx[zcyl_] =          2./ncell[zcyl_];

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double r = c->G.x[r_];
      double z = c->G.x[zcyl_];

      c->S.prim[RHO] = 1.;
      c->S.prim[VV1] = 0.;
      c->S.prim[VV2] = 0.;
      c->S.prim[PPP] = 1.;
      c->S.prim[TR1] = 2.;

      if (z < -0.5 and fabs(r)< 0.2){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.7;
        c->S.prim[PPP] = 1.;
        c->S.prim[TR1] = 1.;

      }
    }

  }

  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  double vIn     = 0.;     // can't be lower than 1 for algo to work
  double vOut    = 0.;     
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n <= ngst; ++n){
      int    iL = n;
      int    iR = ntrack[j]-2-n;
      Itot[j][iL].v = vIn;
      Itot[j][iR].v = vOut;
    }
  }

}


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  // Cell c = Ctot[j][i];

  // double split_dl = 0.05;
  // double merge_dl = 0.0001;
    // careful, dx != dl

  // if (c.G.dx[MV] > split_dl) {
  //   return(split_);
  // }
  // if (c.G.dx[MV] < merge_dl) {
  //   return(merge_);
  // }
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














