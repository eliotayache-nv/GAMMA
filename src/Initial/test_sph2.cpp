/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-05 17:44:03
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 100;
  par->ncell[y_] = 100;
  par->nmax      = 310;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) 20*(i+0.5)/ncell[x_] + 10.;
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
      if (x < 20){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 0.1;
        c->S.cons[NUM_C] = 1.;
      }
      else{
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.0;
        c->S.prim[VV2] = 0.0;
        c->S.prim[PPP] = 0.01;
        c->S.cons[NUM_C] = 2.;
      }
    }

  }

  return 0;

}


int Cell::checkCellForRegrid(){

  // double split_dl = 1;
  // double merge_dl = 0.05;

  // if (G.dx[MV] > split_dl) {
  //   return(split_);
  // }
  // if (G.dx[MV] < merge_dl) {
  //   return(merge_);
  // }
  return(skip_);

}





