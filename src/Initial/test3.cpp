/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-05 14:00:32
*/

#include "../environment.h"
#include "../grid.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 100;
  par->ncell[y_] = 100;
  par->nmax      = 110;    // max number of cells in MV direction
  par->ngst      = 1;

}

int Grid::initialGeometry(){

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) 2.*(i+0.5)/ncell[x_] - 1.;
      c->G.dx[x_] =          2./ncell[x_];
      c->G.x[y_]  = (double) 2.*(j+0.5)/ncell[y_] - 1.0001;
      c->G.dx[y_] =          2./ncell[y_];
      c->G.dV     = c->G.dx[x_]*c->G.dx[y_];
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
      if (y > x and y >= -x){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.7;
        c->S.prim[VV2] = -0.7;
        c->S.prim[PPP] = 1.;
        c->S.cons[NUM_C] = 2.;

      }
      if (y < -x and y > x){
        c->S.prim[RHO] = 0.5;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 1.;
        c->S.cons[NUM_C] = 3.;
      }
      if (y <= x and y < -x){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.7;
        c->S.prim[VV2] = 0.7;
        c->S.prim[PPP] = 1.;
        c->S.cons[NUM_C] = 4.;

      }
      if (y >= -x and y <= x){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 0.01;
        c->S.cons[NUM_C] = 1.;
      }
    }

  }

  return 0;

}


int Cell::checkCellForRegrid(){

  double split_dl = 0.05;
  double merge_dl = 0.002;

  if (G.dx[MV] > split_dl) {
    // printf("split %d %d\n", nde_ind[y_], nde_ind[x_]);
    return(split_);
  }
  if (G.dx[MV] < merge_dl) {
    // printf("merge %d %d\n", nde_ind[y_], nde_ind[x_]);
    return(merge_);
  }
  return(skip_);

}






