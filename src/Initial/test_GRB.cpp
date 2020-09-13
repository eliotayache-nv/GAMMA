/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-13 10:53:57
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

static double Eiso  = 1.e52;   // erg
static double t90   = 100;       // s
static double eta   = 1.e-3;    
static double n0    = 1.e0;    // cm-3
static double lfac0 = 100;
static double theta0= 0.1;     // rad: jet opening angle
static double rmin  = 5.e11;  // cm : begining of the box at startup
static double rmax  = 1.e13;   // cm : end of the box at startup
static double r0    = 4.e12;  // cm : back of shell
static double r1    = 5.e12;   // cm : head of shell
static double dtheta= 0.2;   // rad; grid opening angle

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 300;
  par->ncell[y_] = 150;
  par->nmax      = 320;    // max number of cells in MV direction
  par->ngst      = 2;

}

// ---------------------------------------------------------------------------------------
// GRB-specific functions






// ---------------------------------------------------------------------------------------
// Simulation functions

int Grid::initialGeometry(){
  // log grid
  double logrmin = log(rmin);
  double logrmax = log(rmax);
  double dlogr = (logrmax - logrmin);

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      double r  = (double) dlogr*(i+0.5)/ncell[x_] + logrmin;
      double rL = (double) dlogr*(i  )/ncell[x_] + logrmin;
      double rR = (double) dlogr*(i+1)/ncell[x_] + logrmin;
      c->G.x[x_]  = exp(r);
      c->G.dx[x_] = exp(rR) - exp(rL);
      c->G.x[y_]  = (double) dtheta*(j+0.5)/ncell[y_];
      c->G.dx[y_] = dtheta/ncell[y_];
      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  double rho0 = n0*mp_;
  double p0   = eta*rho0;
  double lfac02 = lfac0*lfac0;
  double vr1  = sqrt(1-1./(lfac0*lfac0));
  double Edot = Eiso/t90;
  double k = GAMMA_/(GAMMA_-1);
  double c2 = c_*c_;

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double r    = c->G.x[r_];
      double th   = c->G.x[t_];
      double r2   = r*r;
      double rho1 = Edot / (4*PI*r2*vr1*c_*lfac02*c2*(1 + eta*(k - 1./lfac02)) );
      double p1   = eta*rho1;

      if (fabs(th) < theta0 and r0 < r and r < r1){
        c->S.prim[RHO] = rho1;
        c->S.prim[VV1] = vr1;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = p1;
        c->S.cons[NUM_C] = 1.;
      }
      else{
        c->S.prim[RHO] = rho0;
        c->S.prim[VV1] = 0.0;
        c->S.prim[VV2] = 0.0;
        c->S.prim[PPP] = p0;
        c->S.cons[NUM_C] = 2.;            
      }
    }
  }

  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  int check_dist = 90;
  double vIn     = 1.5;     // can't be lower than 1 for algo to work
  double vOut    = 1.2;     
  double vInlim  = 0.5;     // threshold to detect back of the ejecta
  double vOutlim = 0.1;

  int moveInner = 1;
  int moveOuter = 0;
  for (int j = 0; j < nde_nax[F1]; ++j){

    Interface Iin  = Itot[j][iLbnd[j]+check_dist];
    Interface Iout = Itot[j][iRbnd[j]-check_dist];

    if (Iin.v  > vInlim)  moveInner = 0;
    if (Iout.v > vOutlim) moveOuter = 1;

  }

  // communicating to all nodes
  double allmoveInner;
  double allmoveOuter;
  MPI_Allreduce(&moveInner, &allmoveInner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&moveOuter, &allmoveOuter, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

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

  // We always want ncells[x_] cells in the y direction
  // We allow dx in [0.1, 10] around target uniform resolution

  Cell c = Ctot[j][i];
  double split_ratio = 10;    // relative to target resolution
  double merge_ratio = 0.1;   // relative to target resolution
  // double r  = c.G.x[MV];
  double dl = c.G.dl[MV];
  double rmin = Ctot[j][iLbnd[j]+1].G.x[r_];
  double rmax = Ctot[j][iRbnd[j]-1].G.x[r_];
  double target_res = (rmax - rmin) / ncell[r_];

  if (dl > split_ratio * target_res) {
    return(split_);
  }
  if (dl < merge_ratio * target_res) {
    return(merge_);
  }
  return(skip_);

}



void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  double rho0 = n0*mp_;
  double p0   = eta*rho0;
  double ratio = 1.e-8;

  *p = fmax(*p, ratio*p0);
  *rho = fmax(*rho, ratio*rho0);

  return;

}


