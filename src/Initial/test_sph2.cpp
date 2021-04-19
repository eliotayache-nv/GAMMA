/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-18 23:04:19
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 20;
  par->ncell[y_] = 20;
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
      if (x < 20 and y<=PI/24.){
        c->S.prim[RHO] = 0.1;
        c->S.prim[VV1] = 0.7;
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


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  // int check_dist = 90;
  double vIn     = 0.;     // can't be lower than 1 for algo to work
  double vOut    = 0.;     
  // double rhoInlim = 1.e3;     // threshold to detect back of the ejecta
  // double vOutlim = 0.1;

  // int moveInner = 1.;
  // int moveOuter = 1.;
  // for (int j = 0; j < nde_nax[F1]; ++j){

  //   Cell      Cin  = Ctot[j][iLbnd[j]+check_dist];
  //   Interface Iout = Itot[j][iRbnd[j]-check_dist];

  //   if (Cin.S.prim[RHO] > rhoInlim)  moveInner = 0;
  //   if (Iout.v > vOutlim) moveOuter = 1;

  // }

  // communicating to all nodes
  double allmoveInner = 1;
  double allmoveOuter = 1;
  // double allmoveInner;
  // double allmoveOuter;
  // MPI_Allreduce(&moveInner, &allmoveInner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  // MPI_Allreduce(&moveOuter, &allmoveOuter, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n <= ngst; ++n){
      int    iL = n;
      int    iR = ntrack[j]-2-n;
      if (allmoveInner) Itot[j][iL].v = vIn;
      if (allmoveOuter) Itot[j][iR].v = vOut;
    }
  }
  // for (int j = 0; j < nde_nax[F1]; ++j){
  //   for (int n = 0; n <= ngst; ++n){
  //     int    iL = n;
  //     int    iR = ntrack[j]-2-n;
  //     Itot[j][iL].v = 0;
  //     Itot[j][iR].v = 0;
  //   }
  // }

}


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell *c = &Ctot[j][i];
  double split_dl = 2;
  double merge_dl = 0.2;

  if (c->G.dx[MV] > split_dl) {
    return(split_);
  }
  if (c->G.dx[MV] < merge_dl) {
    return(merge_);
  }
  return(skip_);

}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value

  // double trac = S.prim[TR1];
  // if (fabs(trac-1 < 0.2)){
  //   *res /= ejecta_regrid_ratio;
  // }
  // if (fabs(trac-2 < 0.5)){
  //   *res /= csm_regrid_ratio;
  // }

}



void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  // double ratio = 1.e-5;

  // *p = fmax(*p, ratio*p0/pNorm);
  // *rho = fmax(*rho, ratio*rho0/rhoNorm);

  return;

}




