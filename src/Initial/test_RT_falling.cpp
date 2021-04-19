/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-11-16 17:48:14
*/

#include "../environment.h"
#include "../grid.h"
#include "../constants.h"
#include <gsl/gsl_roots.h>   // root finding algorithms
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>


double g = 0.1; 
double p0 = 2.5;

void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = 100;
  par->ncell[y_] = 25;
  par->nmax      = 110;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){
  // Careful! When switching moving coordinate, you need to set MV as a function of i
  // and F1 as a functino of j

  double xmin = -0.75;
  double xmax =  0.75;
  double ymin = -0.25;
  double ymax =  0.25;

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      c->G.x[x_]  = (double) (xmax-xmin)*(i+0.5)/ncell[x_] + xmin;
      c->G.dx[x_] =          (xmax-xmin)/ncell[x_];
      c->G.x[y_]  = (double) (ymax-ymin)*(j+0.5)/ncell[y_] + ymin;
      c->G.dx[y_] =          (ymax-ymin)/ncell[y_];

      c->computeAllGeom();
    }
  }
  return 0;

}

int Grid::initialValues(){

  // const gsl_rng_type * T;
  // gsl_rng * r;

  // gsl_rng_env_setup();
  // T = gsl_rng_default;
  // r = gsl_rng_alloc (T);

  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];

      double x = c->G.x[x_];
      double y = c->G.x[y_];

      double rho, trac;
      if (x < 0){
        rho = 1.;
        trac = 1.;
      }
      else {
        rho = 2.;
        trac = 2.;
      }
      c->S.prim[RHO] = rho;
      c->S.prim[TR1] = trac;
      c->S.prim[UU1] = 0.01 * ((1. + cos(4.*PI*y)) * (1. + cos(3.*PI*x))) / 4.;
      c->S.prim[VV2] = 0.;
      c->S.prim[PPP] = p0 - g*rho*x;

    }
  }

  gsl_rng_free (r);
  return 0;

}


void Grid::userKinematics(){

  // // setting lower and higher i boundary interface velocities to zero
  // double vIn     = 0.;     // can't be lower than 1 for algo to work
  // double vOut    = 0.;     
  // for (int j = 0; j < nde_nax[F1]; ++j){
  //   for (int n = 0; n <= ngst; ++n){
  //     int    iL = n;
  //     int    iR = ntrack[j]-2-n;
  //     Itot[j][iL].v = vIn;
  //     Itot[j][iR].v = vOut;
  //   }
  // }

}


void Cell::userSourceTerms(double dt){

  // sources have to be expressed in terms of conserved variables
  double dV  = G.dV;
  double rho = S.prim[RHO];
  double lfac = S.lfac();
  double vx  = S.prim[VV1]/lfac;
  double fg = -lfac*rho*g*dV*dt;
  double Wg = fg*vx;

  S.cons[SS1] += fg;
  S.cons[TAU] += Wg;


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
      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  if (worldrank==worldsize-1){
    for (int j = nde_nax[F1]-ngst; j < nde_nax[F1]; ++j){
      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  // Periodic boundary conditions
  // int rL = 0;             // rank containing left  boundary
  // int rR = worldsize-1;   // rank containing right boundary
  // MPI_Datatype cell_mpi = {0}; 
  // generate_mpi_cell( &cell_mpi );
  // MPI_Barrier(MPI_COMM_WORLD);

  // if (worldrank == 0){
  //   for (int j = 0; j <= jLbnd; ++j){

  //     int jout = j+ngst;
  //     int jin  = j;
  //     int nout = nact[jout]+2*ngst; // number of cells to receive
  //     int nin  = 1;                 // number of cells to send (default is 1)
  //     int tag1 = j;
  //     int tag2 = ngst+j;

  //     MPI_Sendrecv(&nout, 1, MPI_INT, rL, tag1, 
  //                  &nin , 1, MPI_INT, rR, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   

  //     ntrack[jin] = nin;
  //     nact[jin]  = nin-2*ngst;
  //     iRbnd[jin] = nin-ngst; 

  //     s_cell *SCout = new s_cell[nout];
  //     s_cell *SCin  = new s_cell[nin];

  //     for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
  //     MPI_Sendrecv(SCout, nout, cell_mpi, rL, tag2, 
  //                  SCin , nin , cell_mpi, rR, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }          

  //     delete(SCin); delete(SCout);

  //   }
  // }

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  // Cell c = Ctot[j][i];

  // double split_dl = 0.05;
  // double merge_dl = 0.0005;
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














