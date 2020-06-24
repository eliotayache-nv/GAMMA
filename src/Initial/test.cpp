/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-23 16:59:42
*/

#include "../environment.h"
#include "../grid.h"

void loadParams(s_par *par){

  par->tini       = 0.;
  par->ncell[x_] = 10;
  par->ncell[y_] = 10;
  par->ngst      = 1;

}

int Grid :: initialGeometry(){

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double)i/ncell[x_];
      c->G.dl[x_] =        1./ncell[x_];
      c->G.x[y_]  = (double)j/ncell[y_];
      c->G.dl[y_] =        1./ncell[y_];
      c->G.dV     = c->G.dl[x_]*c->G.dl[y_];
    }
  }
  return 0;

}

int Grid :: initialValues(){

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];
      c->S.prim[RHO] = 1;
      c->S.prim[PPP] = 1;
      c->S.prim[UU1] = 0.1;
      c->S.prim[UU2] = 0;
    }
  }
  return 0;

}