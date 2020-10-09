/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-09 09:39:48
*/

// rupert simulation 2: 

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

void loadParams(s_par *par){

  par->tini      = 0.;      // initial time(?)
  par->ncell[x_] = 200;     // number of cells in r direction
  par->ncell[y_] = 10;      // number of cells in theta direction
  par->nmax      = 1000;    // max number of cells in MV direction
  par->ngst      = 2;       // number of ghost cells (?); probably don't change

}

int Grid::initialGeometry(){                              // loop across grid, calculate cell positions and grid spacing, add info to cell struct "c"

// Changes:
// Added settings for changing the grid size (doesn't change number of cells)

  double grid_rsize = 2.;
  double grid_thetasize = 10.*2.*PI/360.;

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) grid_rsize*i/ncell[x_];                  // r position
      c->G.dx[x_] =          grid_rsize/ncell[x_];                    // r grid spacing
      c->G.x[y_]  = (double) grid_thetasize*(j+0.5)/ncell[y_];              // theta position
      c->G.dx[y_] =          grid_thetasize/ncell[y_];                // theta grid spacing

      c->computeAllGeom();
    }
  }
  return 0;
}

int Grid::initialValues(){

  // set up shell in lab frame
  double shl_rcentre_lab = 0.5;      // set r value about which shell is centred (in lab frame)
  double shl_rwidth_lab = 0.2;       // set thickness of shell in lab frame
  
  // set comoving fluid quantities, in the jet and in the shell
  double rho_jet = 0.1;
  double p_jet = 0.1;
  double rho_shl = 1.;
  double p_shl = 1.;

  // set velocities in the lab frame
  double v_shl_lab = 0.5;        // set shell velocity (along r) in lab frame
  double v_jet_lab = 0;           // set jet velocity (along r) in lab frame
  
  // set velocity of the reference frame
  double v_ref = 0.4;             // set velocity of reference frame
  
  // lorentz transformations
  double v_shl_ref, v_jet_ref, lorentz_gamma, shl_rcentre_ref, shl_rwidth_ref;
  lorentz_gamma = 1./sqrt(1-(v_ref*v_ref));                        // calculate lorentz factor
  shl_rcentre_ref = lorentz_gamma*shl_rcentre_lab;              // calculate shell centre position in reference frame
  shl_rwidth_ref = lorentz_gamma*shl_rwidth_lab;                // calculate shell width in reference frame
  v_jet_ref = (v_jet_lab - v_ref)/(1 - v_ref*v_jet_lab);          // calculate jet velocity in reference frame
  v_shl_ref = (v_shl_lab - v_ref)/(1 - v_ref*v_shl_lab);       // calculate shell velocity in reference frame

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[j][i];

      double x = c->G.x[x_];
      double y = c->G.x[y_];

      c->S.prim[RHO] = rho_jet;
      c->S.prim[VV1] = v_jet_ref;
      c->S.prim[VV2] = 0.;
      c->S.prim[PPP] = p_jet;
      c->S.prim[TR1] = 2.;

      if (fabs(x-shl_rcentre_ref) < shl_rwidth_ref/2.){
        c->S.prim[RHO] = rho_shl;
        c->S.prim[VV1] = v_shl_ref;
        c->S.prim[VV2] = 0.;
        c->S.prim[PPP] = p_shl;
        c->S.prim[TR1] = 1.;

      }
    }

  }

  return 0;

}


void Grid::userKinematics(){

  // setting lower and higher i boundary interface velocities to zero
  double vIn     = 0.;     // can't be lower than 1 for algo to work
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

  // Cell c = Ctot[j][i];

  // double split_dl = 0.05;
  // double merge_dl = 0.0001;
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
  UNUSED(*rho);
  UNUSED(*p);

  return;

}














