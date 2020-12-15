/*
* @Author: Eliot Ayache
* @Date:   2020-12-08 21:44:52
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-12-15 10:24:35
*/

/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-08 16:00:45
*/

// rupert simulation 1: 

#include "../../environment.h"
#include "../../grid.h"

double lScale = 100.;        // scaling applied to lengths and velocities (velocities decreased by a factor of 100, so c_sound = 3.5/100 < c_light)

static double rhoScale = 1./(lScale*lScale*lScale);   // rho = M / L3
static double vScale = lScale;                        // v   = L / T
static double aScale = lScale;                        // a   = L / T2
static double pScale = aScale/(lScale*lScale);        // p   = Ma / L2

static double grav        = 0.1/aScale;               // gravity is acceleration
static double grid_xsize  = 1.5/lScale;               // length
static double grid_ysize  = 0.5/lScale;               // length
static double kx          = 3.*M_PI*lScale;           // 1/length
static double ky          = 4.*M_PI*lScale;           // 1/length
static double w0          = 0.01/vScale;              // velocity

static double rho2        = 2./rhoScale;              // density
static double rho1        = 1./rhoScale;              // density
static double P0          = 2.5/pScale;               // pressure

void loadParams(s_par *par){

  par->tini      = 0.;      // initial time(?)
  par->ncell[x_] = 74;     // number of cells in x direction
  par->ncell[y_] = 24;     // number of cells in y direction
  par->nmax      = 300;    // max number of cells in MV direction
  par->ngst      = 2;       // number of ghost cells on each side

}

int Grid::initialGeometry(){                              // loop across grid, calculate cell positions and grid spacing, add info to cell struct "c"

// Changes:
// Added settings for changing the grid size (doesn't change number of cells)

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) grid_xsize*(i+0.5)/ncell[x_] - grid_xsize/2.;   // x position
      c->G.dx[x_] =          grid_xsize/ncell[x_];                // x grid spacing
      c->G.x[y_]  = (double) grid_ysize*(j+0.5)/ncell[y_] - grid_ysize/2.;   // y position
      c->G.dx[y_] =          grid_ysize/ncell[y_];                // y grid spacing

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  if (GAMMA_ != 1.4){
    printf("change GAMMA_ to 1.4\n");
  }
  printf("VI\t= %f\n",VI);

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double x    = c->G.x[x_];
      double y    = c->G.x[y_];

      double vx     = w0*(1+cos(ky*y))*(1+cos(kx*x))/4.;
      // double vx     = 0;

      double rho_t,P_t;

      if (x <= 0){
        rho_t      = rho1;
        P_t        = P0 - grav*rho_t*x;
        c->S.prim[RHO]  = rho_t;
        c->S.prim[VV1]  = vx;
        c->S.prim[VV2]  = 0.;
        c->S.prim[PPP]  = P_t;
        c->S.prim[TR1]  = 1.;
      }
      else{
        rho_t      = rho2;
        P_t        = P0 - grav*rho_t*x;
        c->S.prim[RHO]  = rho_t;
        c->S.prim[VV1]  = vx;
        c->S.prim[VV2]  = 0.;
        c->S.prim[PPP]  = P_t;
        c->S.prim[TR1]  = 2.;
      }

    }
  }
  return 0;

}


void Grid::userKinematics(int it, double t){

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

  double rho = S.prim[RHO];               // get density in cell
  double lfac = S.lfac();                 // calculate lorentz factor
  double D = rho*lfac;                    // calculate lab frame density
  double v_x = S.prim[UU1]/lfac;          // calculate lab frame velocity
  double dV = G.dV;                       // get cell volume; must multiply source terms by this to add to conserved variables

  // calculate source terms; multiply by dt since equation is dU/dt + div.flux = source
  S.cons[SS1] += -D*grav*dV*dt;           // momentum source term
  S.cons[TAU] += -D*v_x*grav*dV*dt;       // energy source term

}

void Grid::userBoundaries(int it, double t){

  // Reflective boundary conditions on top and bottom
  // this should essentially just keep the fluid in the grid.

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ngst; ++i){
      int nt = ntrack[j];
      Ctot[j][i].S.prim[UU1] *= -1;  
      Ctot[j][nt-1-i].S.prim[UU1] *= -1;  
    }
  }

  // keeping outflow BCs on the sides, rather than reflective
  // ideally I'd use periodic BCs, but these are awkward to implement!

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
        Ctot[j][i].G.x[F1] = Ctot[jLbnd+1][iLbnd[jLbnd+1]+1].G.x[F1];
        Ctot[j][i].G.x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] = Itot[jLbnd+1][iLbnd[jLbnd+1]+1].x[F1];
          Itot[j][i].x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1]; 
          Itot[j][i].computedA();
        }
        // if (i==10) printf("%le\n", Itot[j][i].x[F1]);
      }

      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  if (worldrank==worldsize-1){
    for (int j = nde_nax[F1]-ngst; j < nde_nax[F1]; ++j){

      int target_j = nde_nax[F1] - ngst -1 - (j - jRbnd);
      // printf("%d %d %d %d\n", jRbnd, j, target_j, nde_nax[F1]);
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
        Ctot[j][i].G.x[F1] = Ctot[jRbnd-1][iLbnd[jRbnd-1]+1].G.x[F1];
        Ctot[j][i].G.x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] = Itot[jRbnd-1][iLbnd[jRbnd-1]+1].x[F1];
          Itot[j][i].x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
          Itot[j][i].computedA();
        }
        // if (i==10) printf("%le\n", Itot[j][i].x[F1]);
      }
      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
      // exit(20);
  }

  // printf("%le %le\n", Ctot[1][10].S.prim[RHO], Ctot[jRbnd][10].S.prim[RHO]);
  // printf("%le %le\n", Itot[0][10].x[x_], Itot[jRbnd+1][10].x[x_]);

  // UNUSED(it);
  // UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  // Cell c = Ctot[j][i];

  // double split_dl = 0.002;
  // double merge_dl = 0.05;
  //   // careful, dx != dl

  // if (c.G.dx[MV] > split_dl) {
  //   return(split_);
  // }
  // if (c.G.dx[MV] < merge_dl) {
  //   return(merge_);
  // }
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

  // if (*p<1.e-8){ *p = 1.e-8; }
  // if (*rho<1.e-8){ *rho = 1.e-8; }

  return;

}














