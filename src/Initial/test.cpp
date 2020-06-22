/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-18 15:27:24
*/

#include "../environment.h"
#include "../grid.h"

void loadParams(s_par *par){

  par->tini       = 0.;
  par->n_cell[x_] = 10;
  par->n_cell[y_] = 10;
  par->n_gst      = 1;

}

int Grid :: initialGeometry(){

  for (int j = 0; j < n_cell[y_]; ++j){
    for (int i = 0; i < n_cell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double)i/n_cell[x_];
      c->G.dl[x_] =        1./n_cell[x_];
      c->G.x[y_]  = (double)j/n_cell[y_];
      c->G.dl[y_] =        1./n_cell[y_];
      c->G.dV     = c->G.dl[x_]*c->G.dl[y_];
    }
  }
  return 0;

}

int Grid :: initialValues(){

  for (int j = 0; j < n_cell[F1]; ++j){
    for (int i = 0; i < n_cell[MV]; ++i){
      Cell *c = &Cinit[j][i];
      c->S.prim[RHO] = 1;
      c->S.prim[PPP] = 1;
    }
  }
  return 0;

}