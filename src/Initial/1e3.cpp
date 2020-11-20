/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-08 16:00:45
*/

// rupert simulation 1: 

#include "../environment.h"
#include "../grid.h"

static double grid_xsize  = 0.5;
static double grid_ysize  = 0.25;
static double b           = 0.025;
static double a           = 0.05;

static double rhoMove     = 100.;
static double rhoStatic   = 1.;
static double P           = 1.;
static double vx          = 0.05;
static double vy          = 0.05;

void loadParams(s_par *par){

  par->tini      = 0.;      // initial time(?)
  par->ncell[x_] = 200;     // number of cells in x direction
  par->ncell[y_] = 100;     // number of cells in y direction
  par->nmax      = 400;    // max number of cells in MV direction
  par->ngst      = 2;       // number of ghost cells on each side

}

int Grid::initialGeometry(){                              // loop across grid, calculate cell positions and grid spacing, add info to cell struct "c"

// Changes:
// Added settings for changing the grid size (doesn't change number of cells)

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) grid_xsize*i/ncell[x_] - grid_xsize/2.;   // x position
      c->G.dx[x_] =          grid_xsize/ncell[x_];                // x grid spacing
      c->G.x[y_]  = (double) grid_ysize*j/ncell[y_];   // y position
      c->G.dx[y_] =          grid_ysize/ncell[y_];                // y grid spacing

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  printf("Gamma = %f\n",GAMMA_);
  printf("VI = %f\n",VI);

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double x    = c->G.x[x_];
      double y    = c->G.x[y_];

      if ((y < a + b) and (y > a - b) and (x < a + b - grid_xsize/2.) and (x > a - b - grid_xsize/2.)){

        c->S.prim[RHO]  = rhoMove;
        c->S.prim[VV1]  = vx;
        c->S.prim[VV2]  = vy;
        c->S.prim[PPP]  = P;
        c->S.prim[TR1]  = 2.;
      }
      else{
        c->S.prim[RHO]  = rhoStatic;
        c->S.prim[VV1]  = 0.;
        c->S.prim[VV2]  = 0.;
        c->S.prim[PPP]  = P;
        c->S.prim[TR1]  = 1.;
      }

    }
  }
  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  double vIn     = 0.;
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

void Cell::userSourceTerms(double dt){

  UNUSED(dt);

}

void Grid::userBoundaries(int it, double t){

  // Reflective boundary conditions
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ngst; ++i){
      int nt = ntrack[j];
      Ctot[j][i].S.prim[UU1] *= -1;  
      Ctot[j][nt-1-i].S.prim[UU1] *= -1;  
    }
  }

  if (worldrank==0){
    for (int j = 0; j < ngst; ++j){

      int target_j = 2*ngst-1-j;
      ntrack[j] = ntrack[target_j];
      nact[j] = nact[target_j];
      iRbnd[j] = iRbnd[target_j];
      iLbnd[j] = iLbnd[target_j];
      std::copy_n(&Ctot[target_j][0], ntrack[target_j],   &Ctot[j][0]);
      std::copy_n(&Itot[target_j][0], ntrack[target_j]-1, &Itot[j][0]);

      // updating ghost positions, ids and indexes
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1]; 
          Itot[j][i].computedA();
        }
      }

      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  if (worldrank==worldsize-1){
    for (int j = nde_nax[F1]-ngst; j < nde_nax[F1]; ++j){

      int target_j = nde_nax[F1] - ngst -j;
      ntrack[j] = ntrack[target_j];
      nact[j] = nact[target_j];
      iRbnd[j] = iRbnd[target_j];
      iLbnd[j] = iLbnd[target_j];
      std::copy_n(&Ctot[target_j][0], ntrack[target_j],   &Ctot[j][0]);
      std::copy_n(&Itot[target_j][0], ntrack[target_j]-1, &Itot[j][0]);

      // updating ghost positions, ids and indexes
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
          Itot[j][i].computedA();
        }
      }
      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

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
  UNUSED(p);
  UNUSED(rho);

  return;

}














