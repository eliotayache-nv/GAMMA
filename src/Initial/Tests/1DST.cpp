/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-04-17 10:31:28
*/

#include "../../environment.h"
#include "../../grid.h"
#include "../../simu.h"
#include "../../constants.h"

static double xmin = 0.02;
static double xmax = 0.5;
static double xmid = 0.25;
static double Nx = 10000;


void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = Nx;
  par->nmax      = Nx+50;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){

  double x = xmax - xmin;
  for (int i = 0; i < ncell[x_]; ++i){
    Cell *c = &Cinit[i];
    c->G.x[x_]  = (double) x*(i+0.5)/ncell[x_] + xmin;
    c->G.dx[x_] =          x/ncell[x_];
    c->computeAllGeom();
  }
  return 0;

}


int Grid::initialValues(){

  for (int i = 0; i < ncell[MV]; ++i){
    Cell *c = &Cinit[i];

    // spherical
    if (c->G.x[x_] <= xmid){
      c->S.prim[RHO] = 1;
      c->S.prim[VV1] = 0;
      c->S.prim[PPP] = 1;
      c->S.cons[TR1] = 1.;
    }
    if (c->G.x[x_] > xmid){
      c->S.prim[RHO] = 0.1;
      c->S.prim[VV1] = 0.;
      c->S.prim[PPP] = 0.1;
      c->S.cons[TR1] = 2.;
    }

    //cartesian
    // if (c->G.x[x_] <= xmid){
    //   c->S.prim[RHO] = 0.1;
    //   c->S.prim[VV1] = 0.99;
    //   c->S.prim[PPP] = 1.;
    //   c->S.cons[TR1] = 1.;
    // }
    // if (c->G.x[x_] > xmid){
    //   c->S.prim[RHO] = 1;
    //   c->S.prim[VV1] = 0.;
    //   c->S.prim[PPP] = 1.;
    //   c->S.cons[TR1] = 2.;
    // }
  }

  return 0;

}


void Grid::userKinematics(int it, double t){

  UNUSED(it);
  UNUSED(t);

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

  return;

}



void Simu::dataDump(){

  if (it%1000 == 0){ grid.printCols(it, t); }
  // if (it%100 == 0){ grid.printCols(it, t); }

}

void Simu::runInfo(){

  // if ((worldrank == 0) and (it%1 == 0)){ printf("it: %ld time: %le\n", it, t);}
  if ((worldrank == 0) and (it%100 == 0)){ printf("it: %ld time: %le\n", it, t);}

}

void Simu::evalEnd(){

  // if (it > 300){ stop = true; }
  if (t > 0.3){ grid.printCols(it, t); stop = true; } // 3.33e8 BOXFIT simu

}



