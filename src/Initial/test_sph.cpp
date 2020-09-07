/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-07 16:50:55
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 150;
  par->ncell[y_] = 50;
  par->nmax      = 160;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) 20*(i+0.5)/ncell[x_] + 1.;
      c->G.dx[x_] =          20./ncell[x_];
      c->G.x[y_]  = (double) PI/6.*(j+0.5)/ncell[y_];
      c->G.dx[y_] =          PI/6./ncell[y_];
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
      if (x < 3 and fabs(y) < 0.2){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.99;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 0.1;
        c->S.cons[NUM_C] = 1.;
      }
      else{
        c->S.prim[RHO] = 0.01;
        c->S.prim[VV1] = 0.0;
        c->S.prim[VV2] = 0.0;
        c->S.prim[PPP] = 0.01;
        c->S.cons[NUM_C] = 2.;
      }
    }
    Cell *c0 = &Cinit[j][0];
    c0->S.prim[RHO] = 0.1;
    c0->S.prim[VV1] = 0.;
    c0->S.prim[VV2] = 0.;
    c0->S.prim[PPP] = 0.01;
    c0->S.cons[NUM_C] = 1.;
  }

  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n <= ngst; ++n){
      int iL = n;
      int iR = nact[j]-2-n;
      Itot[j][iL].v = 0;
      Itot[j][iR].v = 0;
    }
  }

}


void Grid::userBoundaries(int it){

  UNUSED(it);
  // for (int j = 0; j < nde_nax[F1]; ++j){
  //   for (int i = 0; i <= iLbnd[j]; ++i){
  //     Cell *c = &Ctot[j][i];
  //     double y = c->G.x[y_];
  //     if (fabs(y) < 0.2){
  //       c->S.prim[RHO] = 0.1;
  //       c->S.prim[VV1] = 0.0;
  //       c->S.prim[VV2] = 0.;
  //       c->S.prim[PPP] = 0.01;
  //       c->S.cons[NUM_C] = 1.;
  //     }
  //     else{
  //       c->S.prim[RHO] = 0.1;
  //       c->S.prim[VV1] = 0.0;
  //       c->S.prim[VV2] = 0.0;
  //       c->S.prim[PPP] = 0.01;
  //       c->S.cons[NUM_C] = 2.;
  //     }
  //   }
  // }  

}


int Cell::checkCellForRegrid(){

  double split_dl = 0.3;
  double merge_dl = 0.01;

  if (G.dx[MV] > split_dl) {
    return(split_);
  }
  if (G.dx[MV] < merge_dl) {
    return(merge_);
  }
  return(skip_);

}





