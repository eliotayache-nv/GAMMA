/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-24 18:49:31
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

static double Eiso  = 1.e51;   // erg
static double t90   = 10;     // s
static double eta   = 1.e-3;    
static double n0    = 1.e0;    // cm-3
static double lfac0 = 100;
static double theta0= 0.1;     // rad: jet opening angle
static double r0    = 3.e11;   // cm : back of shell
static double r1    = 6.e11;   // cm : head of shell

static double zmin  = 2.e11;   // cm : begining of the box at startup
static double zmax  = 7.e11;   // cm : end of the box at startup
static double rmin  = 0.;   // cm : begining of the box at startup
static double rmax  = 5 *zmax*tan(theta0);   // cm : end of the box at startup

static double rho0  = n0*mp_;
static double p0    = rho0*c_*c_*eta;

static double rhoNorm = rho0;
static double lNorm = c_;
static double vNorm = c_;
static double pNorm = rhoNorm*vNorm*vNorm;

static double ejecta_regrid_ratio = 0.1;
static double csm_regrid_ratio    = 1;

void loadParams(s_par *par){

  par->tini         = 0.;
  par->ncell[zcyl_] = 100;
  par->ncell[rcyl_] = 100;
  par->nmax         = 1000;    // max number of cells in MV direction
  par->ngst         = 2;

}

// ---------------------------------------------------------------------------------------
// GRB-specific functions




// ---------------------------------------------------------------------------------------
// Simulation functions

int Grid::initialGeometry(){

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      double delta_z = zmax - zmin;
      double dz      = delta_z / ncell[zcyl_];
      double z       = zmin + (i+0.5) * dz;
      double delta_r = rmax - rmin;
      double dr      = delta_r / ncell[rcyl_];
      double r       = rmin + (j+0.5) * dr;      
      Cell *c = &Cinit[j][i];

      c->G.x[rcyl_]  = r / lNorm;
      c->G.dx[rcyl_] = dr / lNorm;
      c->G.x[zcyl_]  = z / lNorm;;
      c->G.dx[zcyl_] = dz / lNorm;;
      c->computeAllGeom();

    }
  }
  return 0;

}

int Grid::initialValues(){

  double lfac02 = lfac0*lfac0;
  double vr1  = sqrt(1-1./(lfac02)) * c_;
  double Edot = Eiso/t90;
  double k = GAMMA_/(GAMMA_-1);
  double c2 = c_*c_;

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double rcyl = c->G.x[rcyl_]*lNorm; // deNormalise for calculation
      double zcyl = c->G.x[zcyl_]*lNorm;
      double r2 = rcyl*rcyl + zcyl*zcyl;  // r2 in spherical
      double r = sqrt(r2);                 // r  in spherical
      double th = atan(rcyl/zcyl);

      double rho1 = Edot / (4*PI*r2*vr1*lfac02*c2*(1 + eta*(k - 1./lfac02)) );
      double p1   = eta*rho1*c2;

      // normalisation
      double rho0_norm = rho0/rhoNorm;
      double p0_norm = p0/(pNorm);
      double rho1_norm = rho1/rhoNorm;
      double p1_norm = p1/(pNorm);
      double vr1_norm = vr1/vNorm;

      if (fabs(th) < theta0 and r0 < r and r < r1){
        c->S.prim[RHO] = rho1_norm;
        c->S.prim[VV1] = vr1_norm * zcyl / r;
        c->S.prim[VV2] = vr1_norm * rcyl / r;
        c->S.prim[PPP] = p1_norm;
        c->S.prim[TR1] = 1.;
      }
      else{
        c->S.prim[RHO] = rho0_norm;
        c->S.prim[VV1] = 0.0;
        c->S.prim[VV2] = 0.0;
        c->S.prim[PPP] = p0_norm;
        c->S.prim[TR1] = 2.;            
      }
    }
  }

  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  // int check_dist = 90;
  double vIn     = 0.8;     // can't be lower than 1 for algo to work
  double vOut    = 1.2;     
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


void Cell::userSourceTerms(double dt){

  // // sources have to be expressed in terms of conserved variables
  // double y   = G.x[y_]; 
  // double dV  = G.dV;
  // double rho = S.prim[RHO];
  // double fg = -rho*g*dV*dt;

  // S.cons[SS1] += fg;

}


int Grid::checkCellForRegrid(int j, int i){

  // We always want ncells[x_] cells in the y direction
  // We allow dx in [0.1, 10] around target uniform resolution

  Cell c = Ctot[j][i];

  double trac = c.S.prim[TR1];
  double split_ratio = 5;     // relative to target aspect ratio
  double merge_ratio = 0.1;    // relative to target aspect ratio
  double target_ar = 1;
  double dr  = c.G.dx[rcyl_];
  double dz  = c.G.dx[zcyl_];
  double ar  = dz / dr;  // aspect ratio

  if (fabs(trac-1) < 0.2){
    split_ratio *= ejecta_regrid_ratio;   // smaller cells in ejecta
    merge_ratio *= ejecta_regrid_ratio;
  }

  if (fabs(trac-2) < 0.5){
    split_ratio *= csm_regrid_ratio;   // bigger cells in csm
    merge_ratio *= csm_regrid_ratio;
  }

  if (ar > split_ratio * target_ar) {
    // printf("%d %d\n", j, i);
    return(split_);
  }
  if (ar < merge_ratio * target_ar) {
    return(merge_);
  }
  return(skip_);

}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value

  double trac = S.prim[TR1];
  if (fabs(trac-1 < 0.2)){
    *res /= ejecta_regrid_ratio;
  }
  if (fabs(trac-2 < 0.5)){
    *res /= csm_regrid_ratio;
  }

}



void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  double ratio = 1.e-5;

  *p = fmax(*p, ratio*p0/pNorm);
  *rho = fmax(*rho, ratio*rho0/rhoNorm);

  return;

}


