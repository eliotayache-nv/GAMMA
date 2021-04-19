/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-08 16:00:45
*/

// rupert simulation 1: 

#include "../environment.h"
#include "../grid.h"

void loadParams(s_par *par){

  par->tini      = 0.;      // initial time(?)
  par->ncell[x_] = 200;     // number of cells in x direction
  par->ncell[y_] = 10;     // number of cells in y direction
  par->nmax      = 1000;    // max number of cells in MV direction
  par->ngst      = 2;       // number of ghost cells (?); probably don't change

}

int Grid::initialGeometry(){                              // loop across grid, calculate cell positions and grid spacing, add info to cell struct "c"

// Changes:
// Added settings for changing the grid size (doesn't change number of cells)

  double grid_xsize = 2.;
  double grid_ysize = 1.;

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) grid_xsize*i/ncell[x_] - grid_xsize/2.;   // x position
      c->G.dx[x_] =          grid_xsize/ncell[x_];                // x grid spacing
      c->G.x[y_]  = (double) grid_ysize*j/ncell[y_] - grid_xsize/2.;   // y position
      c->G.dx[y_] =          grid_ysize/ncell[y_];                // y grid spacing

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  // set up slab in lab frame
  double slab_xcentre_lab = 0.;       // set x value about which slab is centred (in lab frame)
  double slab_xwidth_lab = 0.2;   // set thickness of slab in lab frame
  
  // set comoving fluid quantities, in the environment and in the slab
  double rho_environment = 0.1;
  double p_environment = 0.1;
  double rho_slab = 1.;           // set comoving quantities in the slab
  double p_slab = 1.;

  // set velocities in the lab frame
  double v_slab_lab = -0.9;        // set slab velocity (along x) in lab frame
  double v_env_lab = 0;           // set environment velocity (along x) in lab frame
  
  // set velocity of the reference frame
  double v_ref = -0.9;             // set velocity of reference frame
  
  // lorentz transformations
  double v_slab_ref, v_env_ref, lorentz_gamma, slab_xcentre_ref, slab_xwidth_ref;
  lorentz_gamma = 1/sqrt(1-(v_ref*v_ref));                        // calculate lorentz factor
  slab_xcentre_ref = lorentz_gamma*slab_xcentre_lab;              // calculate slab centre position in reference frame
  slab_xwidth_ref = lorentz_gamma*slab_xwidth_lab;                // calculate slab width in reference frame
  v_env_ref = (v_env_lab - v_ref)/(1 - v_ref*v_env_lab);          // calculate environment velocity in reference frame
  v_slab_ref = (v_slab_lab - v_ref)/(1 - v_ref*v_slab_lab);       // calculate slab velocity in reference frame

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double x = c->G.x[x_];
      double y = c->G.x[y_];

      c->S.prim[RHO] = rho_environment;
      c->S.prim[VV1] = v_env_ref;
      c->S.prim[VV2] = 0.;
      c->S.prim[PPP] = p_environment;
      c->S.prim[TR1] = 2.;

      if (fabs(x-slab_xcentre_ref) < slab_xwidth_ref/2.){
        c->S.prim[RHO] = rho_slab;
        c->S.prim[VV1] = v_slab_ref;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = p_slab;
        c->S.prim[TR1] = 1.;

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


void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell c = Ctot[j][i];

  double split_dl = 0.05;
  double merge_dl = 0.0001;
    // careful, dx != dl

  if (c.G.dx[MV] > split_dl) {
    return(split_);
  }
  if (c.G.dx[MV] < merge_dl) {
    return(merge_);
  }
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














