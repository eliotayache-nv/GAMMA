/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-11 19:10:00
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

static double Eiso  = 1.e50;   // erg
static double t90   = 2;       // s
static double eta   = 1.e-3;    
static double n0    = 1.e0;    // cm-3
static double lfac0 = 100;
static double theta0= 0.1;     // rad: jet opening angle
static double rmin  = 5.e9;  // cm : begining of the box at startup
static double rmax  = 6.e11;   // cm : end of the box at startup
static double r0    = 4.e10;  // cm : back of shell
static double r1    = 10.e10;   // cm : head of shell
static double dtheta= (PI/6.);   // rad; grid opening angle

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 400;
  par->ncell[y_] = 100;
  par->nmax      = 420;    // max number of cells in MV direction
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
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n <= ngst; ++n){
      int iL = n;
      int iR = nact[j]-2-n;
      Itot[j][iL].v = 0;
      Itot[j][iR].v = 0;
    }
  }

}


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);
  // double rho0 = n0*mp_;
  // double p0   = eta*rho0;
  // double lfac02 = lfac0*lfac0;
  // double vr1  = sqrt(1-1./(lfac0*lfac0));
  // double Edot = Eiso/t90;
  // double k = GAMMA_/(GAMMA_-1);
  // double c2 = c_*c_;

  // for (int j = 0; j < nde_nax[F1]; ++j){
  //   for (int i = 0; i <= iLbnd[j]; ++i){
  //     Cell *c = &Ctot[j][i];

  //     double r    = c->G.x[r_];
  //     double th   = c->G.x[t_];
  //     double r2   = r*r;
  //     double rho1 = Edot / (4*PI*r2*vr1*c_*lfac02*c2*(1 + eta*(k - 1./lfac02)) );
  //     double p1   = eta*rho1;

  //     if (fabs(th) < theta0 and t < t90*c_){
  //       c->S.prim[RHO] = rho1;
  //       c->S.prim[VV1] = vr1;
  //       c->S.prim[VV2] = 0.;
  //       c->S.prim[PPP] = p1;
  //       c->S.cons[NUM_C] = 1.;
  //     }
  //     else if (fabs(th) < theta0 and t < 10*t90*c_) {
  //       c->S.prim[RHO] = rho0;
  //       c->S.prim[VV1] = vr1;
  //       c->S.prim[VV2] = 0.0;
  //       c->S.prim[PPP] = p0;
  //       c->S.cons[NUM_C] = 2.;    
  //     }
  //     else if (fabs(th) > theta0) {
  //       c->S.prim[RHO] = rho0;
  //       c->S.prim[VV1] = 0.0;
  //       c->S.prim[VV2] = 0.0;
  //       c->S.prim[PPP] = p0;
  //       c->S.cons[NUM_C] = 2.;    
  //     }
  //   }
  // }  

}


int Cell::checkCellForRegrid(){

  double split_ratio = 1.e-1;
  double merge_ratio = 1.e-3;
  double r  = G.x[MV];
  double dl = G.dl[MV];
  double split_dl = split_ratio * r;
  double merge_dl = merge_ratio * r;

  if (dl > split_dl) {
    return(split_);
  }
  if (dl < merge_dl) {
    return(merge_);
  }
  return(skip_);

}



void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  double rho0 = n0*mp_;
  double p0   = eta*rho0;
  double ratio = 1.e-5;

  *p = fmax(*p, ratio*p0);
  *rho = fmax(*rho, ratio*rho0);

  return;

}


