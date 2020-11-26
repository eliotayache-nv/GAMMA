/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-11-19 11:14:33
*/

#include "../environment.h"
#include "../grid.h"
#include "../cell.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 100;
  par->ncell[y_] = 1;
  par->nmax      = 110;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){

  double xmin = -1;
  double xmax = 1;
  double ymin = 0;
  double ymax = 1;
  double x = xmax - xmin;
  double y = ymax - ymin;
  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) x*(i+0.5)/ncell[x_] + xmin;
      c->G.dx[x_] =          x/ncell[x_];
      c->G.x[y_]  = (double) y*(j+0.5)/ncell[y_] + ymin;
      c->G.dx[y_] =          y/ncell[y_];
      c->G.dV     = c->G.dx[x_]*c->G.dx[y_];
    }
  }
  return 0;

}

int Grid::initialValues(){

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      if (c->G.x[x_] >= 0){
        c->S.prim[RHO] = 1.;
        c->S.prim[VV1] = 0.5;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 0.1;
        c->S.cons[NUM_C] = 1.;
      }
      if (c->G.x[x_] < 0){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = 1.;
        c->S.cons[NUM_C] = 2.;
      }
    }

  }

  return 0;

}


void Grid::userKinematics(){

}

void Cell::userSourceTerms(double dt){

  UNUSED(dt);

}

void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  UNUSED(j);
  UNUSED(i);

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
  UNUSED(*p);

  if (*rho < 0.05) *rho = 0.05;

  return;

}


