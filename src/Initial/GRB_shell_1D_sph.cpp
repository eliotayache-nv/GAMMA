/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-14 10:13:49
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"


static double E     = 1.e53;      // isotropic equivalent energy injection, erg
static double t0    = 0.;         // explosion start time, s
static double t1    = 100.;       // explosion end time, s
static double t_ini = 200.;       // start time of dynamic simulation, s
static double lfac  = 100.;       // Lorentz factor of shell
static double eta   = 1.e-3;      // eta = p/(rho*c^2)
static double n0    = 1.;         // CMB number density, cm-3

static double rmin  = 1.e12;      // box lower boundary at t_ini, cm
static double rmax  = 1.e13;      // box upper boundary at t_ini, cm
static double theta = 0.2;        // simulation angle, rad

static double rho0  = n0*mp_;     // CBM mass density, g.cm-3

static double Edot  = E/(t1-t0);  // energy release p.u. time
static double lfac2 = lfac*lfac;  // Lorentz factor squared

// normalisation constants:
static double rhoNorm = rho0;                 // density normalised to CBM density
static double lNorm = c_;                     // distance normalised to c
static double vNorm = c_;                     // velocity normalised to c
static double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

// calculate R0 and R1
static double r0    = (sqrt(1-(1/lfac2)))*c_*(t_ini-t1-t0);   // initial shell inner radius
static double r1    = (sqrt(1-(1/lfac2)))*c_*(t_ini-t0);      // initial shell outer radius

void loadParams(s_par *par){

  par->tini      = t_ini;         // initial time
  par->ncell[x_] = 1000;          // number of cells in r direction
  par->ncell[y_] = 10;             // number of cells in theta direction
  par->nmax      = 1010;          // max number of cells in MV direction
  par->ngst      = 2;             // number of ghost cells (?); probably don't change

}

int Grid::initialGeometry(){                              // loop across grid, calculate cell positions and grid spacing, add info to cell struct "c"

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      
      double r  = (double) (rmax-rmin)*(i+0.5)/ncell[x_] + rmin;
      double dr = (rmax-rmin)/ncell[x_];
      c->G.x[x_]  = r/lNorm;
      c->G.dx[x_] = dr/lNorm;
      c->G.x[y_]  = (j+0.5)*theta/ncell[y_];
      c->G.dx[y_] = theta/ncell[y_];

      c->computeAllGeom();
    }
  }
  return 0;
}

int Grid::initialValues(){

  printf("r0 = %e\n",r0/lNorm);
  printf("r1 = %e\n",r1/lNorm);
  double k        = GAMMA_/(GAMMA_-1);
  double c2       = c_*c_;
  double p0       = eta*rho0*c2;

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double r = c->G.x[x_] * lNorm;    // denormalise for density calculation
      double th = c->G.x[y_];
      double r2 = r*r;
      //printf("r0 = %e\t",r0);
      //printf("r = %e\t",r);
      //printf("r1 = %e\t",r1);
      double v1_r = sqrt(1.-(1./lfac2))*c_;             // calculate shell outflow velocity at r
      double rho1_r = Edot / (4.*PI*r2*v1_r*lfac2*c2*(1. + eta*(k - 1./lfac2)) );  // calculate shell density at r
      double p1_r = eta*rho1_r*c2;                  // calculate pressure at r

      if (r0 < r and r < r1){
        c->S.prim[RHO] = rho1_r/rhoNorm;
        c->S.prim[VV1] = v1_r/vNorm;
        c->S.prim[VV2] = 0;
        c->S.prim[PPP] = p1_r/pNorm;
        c->S.prim[TR1] = 2.;
      }

      else{
        c->S.prim[RHO] = rho0/rhoNorm;
        c->S.prim[VV1] = 0;
        c->S.prim[VV2] = 0;
        c->S.prim[PPP] = p0/pNorm;
        c->S.prim[TR1] = 1.;
      }
      //printf("rho = %e\t",c->S.prim[RHO]);
      //printf("p = %e\t",c->S.prim[PPP]);
      //printf("tracer = %e\n",c->S.prim[TR1]);
    }

  }

  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities
  double vIn     = 0.;
  double vOut    = 0.5;     
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n <= ngst; ++n){
      int    iL = n;
      int    iR = ntrack[j]-2-n;
      Itot[j][iL].v = vIn;
      Itot[j][iR].v = vOut;
    }
  }

}


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  //Cell c = Ctot[j][i];
//
  ////double split_dl = 0.05;
  ////double merge_dl = 0.0001;
  //  // careful, dx != dl
//
  //if (c.G.dx[MV] > split_dl) {
  //  return(split_);
  //}
  //if (c.G.dx[MV] < merge_dl) {
  //  return(merge_);
  //}
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
  UNUSED(*rho);
  UNUSED(*p);

  return;

}














