/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-10 17:55:37
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

static double Eiso  = 1.e51;   // erg
static double t90   = 2;       // s
static double eta   = 1.e-3;    
static double n0    = 1.e0;    // cm-3
static double lfac0 = 100;
static double theta0= 0.1;     // rad: jet opening angle
static double r0    = 1.e12;   // cm : begining of the box at startup
static double r1    = 1.e13;   // cm : end of the box at startup
static double dtheta= (PI/6.);   // rad; grid opening angle

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 300;
  par->ncell[y_] = 100;
  par->nmax      = 320;    // max number of cells in MV direction
  par->ngst      = 2;

}

// ---------------------------------------------------------------------------------------
// GRB-specific functions






// ---------------------------------------------------------------------------------------
// Simulation functions

int Grid::initialGeometry(){

  double dr = r1 - r0;
  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) dr*(i+0.5)/ncell[x_] + r0;
      c->G.dx[x_] =          dr/ncell[x_];
      c->G.x[y_]  = (double) dtheta*(j+0.5)/ncell[y_];
      c->G.dx[y_] =          dtheta/ncell[y_];
      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  double rho0 = n0*mp_;
  double p0   = eta*rho0 / (c_*c_);
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
      double rho1 = Edot / (4*PI*r2*vr1*lfac02*c2*(1 + eta*(k - 1./lfac02)) );
      double p1   = eta*rho1;

      if (fabs(th) < theta0 and i==0){
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
  double rho0 = n0*mp_;
  double p0   = eta*rho0 / (c_*c_); // to normalise
  double lfac02 = lfac0*lfac0;
  double vr1  = sqrt(1-1./(lfac0*lfac0));
  double Edot = Eiso/t90;
  double k = GAMMA_/(GAMMA_-1);
  double c2 = c_*c_;

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i <= iLbnd[j]; ++i){
      Cell *c = &Ctot[j][i];

      double r    = c->G.x[r_];
      double th   = c->G.x[t_];
      double r2   = r*r;
      double rho1 = Edot / (4*PI*r2*vr1*lfac02*c2*(1 + eta*(k - 1./lfac02)) );
      double p1   = eta*rho1;

      if (fabs(th) < theta0 and t<t90){
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





