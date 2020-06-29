/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-29 10:21:21
*/

#include "../environment.h"
#include "../grid.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 40;
  par->ncell[y_] = 20;
  par->nmax      = 80;    // max number of cells in MV direction
  par->ngst      = 1;

}

int Grid::initialGeometry(){

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) 2.*(i+0.5)/ncell[x_] - 1.;
      c->G.dl[x_] =          2./ncell[x_];
      c->G.x[y_]  = (double) 2.*(j+0.5)/ncell[y_] - 1.;
      c->G.dl[y_] =          2./ncell[y_];
      c->G.dV     = c->G.dl[x_]*c->G.dl[y_];
    }
  }
  return 0;

}

int Grid::initialValues(){

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      if (c->G.x[x_] >= 0 and c->G.x[y_] >= 0){
        c->S.prim[RHO] = 0.1;
        c->S.prim[UU1] = 0.;
        c->S.prim[UU2] = 0.;
        c->S.prim[PPP] = 0.01;
      }
      if (c->G.x[x_] < 0 and c->G.x[y_] >= 0){
        c->S.prim[RHO] = 0.1;
        c->S.prim[UU1] = 0.99;
        c->S.prim[UU2] = 0.;
        c->S.prim[PPP] = 1.;
      }
      if (c->G.x[x_] < 0 and c->G.x[y_] < 0){
        c->S.prim[RHO] = 0.5;
        c->S.prim[UU1] = 0.;
        c->S.prim[UU2] = 0.;
        c->S.prim[PPP] = 1.;
      }
      if (c->G.x[x_] >= 0 and c->G.x[y_] < 0){
        c->S.prim[RHO] = 0.1;
        c->S.prim[UU1] = 0.;
        c->S.prim[UU2] = 0.99;
        c->S.prim[PPP] = 1.;
      }
    }
  }

  return 0;

}