/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-04-11 17:23:33
*/

#include "../environment.h"
#include "../grid.h"
#include "../interface.h"
#include "../cell.h"
#include "../array_tools.h"
#include "../mpisetup.h"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <cstring>


template <> Interface* array_1d_nogst<Interface>(Interface* in, const int ngst){

  Interface* p;
  p = &in[ngst]; 
  return p;

}

// static int splitGrid(int ncell, int ngst, int *origin){
//   /*
//   int ncell:   total number of cells (non including ghost boundaries)
//   int ngst:    size of ghost regions.
//   int *origin:  coordinates of C[0][0] in the simulation domain
//   */
//   int base = ncell/worldsize;
//   int rest = ncell%worldsize;
//   int nde_nax = base;
//   if (worldrank<rest) nde_nax++;
//   nde_nax+=2*ngst;

//   origin[MV] = 0;
//   origin[F1] = worldrank*base;
//   if (worldrank <rest) origin[F1]+=worldrank;
//   else origin[F1]+=rest;

//   return nde_nax;

// }

void Grid::initialise(s_par par){

  for (int d = 0; d < NUM_D; ++d){
    ncell[d] = par.ncell[d];
    ngst     = par.ngst;
    nax[d]   = ncell[d]+2*ngst;
  }
  nax[MV]    = par.nmax+2*ngst;  // overriding with max number of cells
  nsimu = 1; for (int d = 0; d < NUM_D; ++d){ nsimu *= ncell[d]; }
  nde_nax[MV]   = nax[MV];
  origin[MV]    = 0;
  nde_ncell[MV] = ncell[MV];  // max number of active cells in mov dim
  nde_ntot      = nde_nax[MV];

  // active cells and ghost cells
  nact = nde_ncell[MV];
  ntrack = nact+2*ngst;
  iLbnd = ngst-1;
  iRbnd = nde_ncell[MV]+ngst;
  ibig   = 0;
  ismall = 0;
  
  Cinit = array_1d<Cell>(nax[MV]-2*ngst); // Cinit need to have the empty cells...
  Ctot  = array_1d<Cell>(nde_nax[MV]);                       // ...for restart
  Itot  = array_1d<Interface>(nde_nax[MV]-1);
  I     = array_1d_nogst<Interface>(Itot, ngst); // removes ghost to ghost
  C     = array_1d_nogst<Cell>(Ctot, ngst);     

}


void Grid::assignId(int ind[NUM_D]){

  Cell *c = &Ctot[ind[0]];
  c->nde_id = ind[0];
  for (int d = 0; d < NUM_D; ++d){ c->nde_ind[d] = ind[d]; }

}


// void Grid::mpi_exchangeGhostTracks(){

//   int rL = worldrank-1;
//   int rR = worldrank+1;
//   MPI_Datatype cell_mpi = {0}; 
//   generate_mpi_cell( &cell_mpi );
//   MPI_Barrier(MPI_COMM_WORLD);

//   // for interfaces I just need to exchange the positions. That can be done in the cells
//   // directly or I can just exchange a vector of positions

//   // sending to lower node, receiving from higher node
//   for (int j = 0; j < ngst; ++j){
//     int jout = j+ngst;
//     int jin  = jRbnd+j;
//     int nout = nact[jout]+2*ngst; // number of cells to receive
//     int nin  = 1;                 // number of cells to send (default is 1)
//     int tag1 = j;
//     int tag2 = ngst+j;

//     if (worldrank==0){
//       MPI_Recv(&nin , 1, MPI_INT, rR, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
//     } else if (worldrank==worldsize-1){
//       MPI_Send(&nout, 1, MPI_INT, rL, tag1, MPI_COMM_WORLD);      
//     } else {
//       MPI_Sendrecv(&nout, 1, MPI_INT, rL, tag1, 
//                    &nin , 1, MPI_INT, rR, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
//     }

//     // updating ntrack, nact, iRbnd (iLbnd doesn't need to be updated)
//     if (worldrank!=worldsize-1){
//       ntrack[jin] = nin;
//       nact[jin]  = nin-2*ngst;
//       iRbnd[jin] = nin-ngst; 
//     }

//     s_cell *SCout = new s_cell[nout];
//     s_cell *SCin  = new s_cell[nin];
//     if (worldrank==0){
//       MPI_Recv(SCin , nin , cell_mpi, rR, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }    
//     } 
//     else if (worldrank==worldsize-1){
//       for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
//       MPI_Send(SCout, nout, cell_mpi, rL, tag2, MPI_COMM_WORLD);      
//     } 
//     else {
//       for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
//       MPI_Sendrecv(SCout, nout, cell_mpi, rL, tag2, 
//                    SCin , nin , cell_mpi, rR, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }          
//     }
//     delete(SCin); delete(SCout);
//   }

//   // sending to higher node, receiving from lower node
//   for (int j = 0; j < ngst; ++j){
//     int jout = jRbnd-ngst+j;
//     int jin  = j;
//     int nout = nact[jout]+2*ngst; // number of cells to receive
//     int nin  = 1;                 // number of cells to send (default is 1)
//     int tag1 = j;
//     int tag2 = ngst+j;

//     if (worldrank==worldsize-1){
//       MPI_Recv(&nin , 1, MPI_INT, rL, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
//     } else if (worldrank==0){
//       MPI_Send(&nout, 1, MPI_INT, rR, tag1, MPI_COMM_WORLD);      
//     } else {
//       MPI_Sendrecv(&nout, 1, MPI_INT, rR, tag1, 
//                    &nin , 1, MPI_INT, rL, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }

//     // updating ntrack, nact, iRbnd (iLbnd doesn't need to be updated)
//     if (worldrank!=0){
//       ntrack[jin] = nin;
//       nact[jin]  = nin-2*ngst;
//       iRbnd[jin] = nin-ngst;
//     }

//     s_cell *SCout = new s_cell[nout];
//     s_cell *SCin  = new s_cell[nin];
//     if (worldrank==worldsize-1){
//       MPI_Recv(SCin , nin , cell_mpi, rL, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }
//     } 
//     else if (worldrank==0){
//       for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
//       MPI_Send(SCout, nout, cell_mpi, rR, tag2, MPI_COMM_WORLD);      
//     } 
//     else {
//       for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
//       MPI_Sendrecv(SCout, nout, cell_mpi, rR, tag2, 
//                    SCin , nin , cell_mpi, rL, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }          
//     }
//     delete(SCin); delete(SCout);
//   }

//   MPI_Barrier(MPI_COMM_WORLD);

//   for (int j = 0; j < ngst; ++j){ // looping over ghost tracks
//     for (int i = 0; i < nde_nax[MV]; ++i){

//       // updating interface positions from communicated cells
//       if (i != nde_nax[MV]-1){ 
//         int jn = nde_nax[F1]-1-j;
//         Itot[j][i].x[MV]  = Ctot[j][i].G.x[MV] + Ctot[j][i].G.dx[MV]/2.;
//         Itot[jn][i].x[MV] = Ctot[jn][i].G.x[MV] + Ctot[jn][i].G.dx[MV]/2.;
//         Itot[j][i].x[F1]  = Ctot[j][i].G.x[F1];
//         Itot[jn][i].x[F1] = Ctot[jn][i].G.x[F1];
//       }

//       // updating IDs
//       int indL[] = {j,i};
//       int indR[] = {nde_nax[F1]-1-j,i};
//       assignId(indL);
//       assignId(indR);
//     }
//   }

// }


void Grid::updateGhosts(int it, double t){

  for (int i = 0; i <= iLbnd; ++i){
    Ctot[i] = Ctot[iLbnd+1];
    Ctot[i].G.x[MV] -= (iLbnd-i+1) * Ctot[iLbnd+1].G.dx[MV];
    Ctot[i].computeAllGeom();
    int ind[] = {i};
    assignId(ind);
  }
  for (int i = 0; i <= iLbnd-1; ++i){
    Itot[i] = Itot[iLbnd+1];
    Itot[i].x[MV] -= (iLbnd-i+1) * Ctot[iLbnd+1].G.dx[MV];
    Itot[i].computedA();
  }

  for (int i = iRbnd; i < ntrack; ++i){
    Ctot[i] = Ctot[iRbnd-1];
    Ctot[i].G.x[MV] += (i-iRbnd+1) * Ctot[iRbnd-1].G.dx[MV];
    Ctot[i].computeAllGeom();
    int ind[] = {i};
    assignId(ind);
  }
  for (int i = iRbnd; i < ntrack-1; ++i){
    Itot[i] = Itot[iRbnd-1];
    Itot[i].x[MV] += (i-iRbnd) * Ctot[iRbnd-1].G.dx[MV];
    Itot[i].computedA();
  }
  userBoundaries(it, t); // overiding with user-specific boundary conditions

}

#if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_

  static double minmod(double a, double b){

    // double theta = 1.5;
    // double r = a/b;
    if (fabs(a) < fabs(b) && a*b > 0){ return(a); }
    if (fabs(b) < fabs(a) && a*b > 0){ return(b); }
    // double lim = fmax(0, fmin(theta*r, fmin((1+r)/2., theta)));
    // return(b*lim);
    return(0);

  }

  static void grad(Cell cL, Cell cR, int dim, double *grad){

    double xL = cL.G.cen[dim];
    double xR = cR.G.cen[dim];

    for (int q = 0; q < NUM_Q; ++q){
      // reconstruction using primitive variables
      double qL = cL.S.prim[q];
      double qR = cR.S.prim[q];
      grad[q] = (qR - qL) / (xR-xL);
    }

  }

  static void gradcons(Cell cL, Cell cR, int dim, double *grad){

    double xL = cL.G.cen[dim];
    double xR = cR.G.cen[dim];

    for (int q = 0; q < NUM_Q; ++q){
      // conserved variables when splitting zones
      double qL = cL.S.cons[q];
      double qR = cR.S.cons[q];
      grad[q] = (qR - qL) / (xR-xL);
    }

  }

  // static Cell* findNeighbor(Cell *c, int side, double x_mv, Grid *g){

  //   Cell *c_out;
  //   Cell **Ctot = g->Ctot;
  //   int n = 0;
  //   int idn;
  //   int n_neigh = c->neigh[F1][side].size();
  //   double xn_mv, dln_mv;
  //   do {
  //     idn = c->neigh[F1][side][n];
  //     c_out  = &Ctot[0][idn];
  //     xn_mv  = c_out->G.x[MV];
  //     dln_mv = c_out->G.dx[MV];
  //     n++;
  //   } while (xn_mv + 0.5*dln_mv < x_mv and n < n_neigh);
  //     // if no aligned neighbor, we return state from closest left neighbor to the left

  //   return c_out;

  // void Grid::computeTransGradients(int j, int i){

  //   Cell *c0 = &Ctot[j][i];
  //   double x0 = c0->G.x[MV];
  //   double xL0 = x0 - c0->G.dx[MV]/2.;
  //   double xR0 = x0 + c0->G.dx[MV]/2.;
  //   const int n_neighs = (int) (c0->neigh[F1][0].size() + c0->neigh[F1][1].size());
  //   double grads[n_neighs][NUM_Q];
  //   double grad_avg[NUM_Q] = {0};
  //   double dAtot = 0;

  //   int ind = 0;
  //   for (int d = 0; d < 2; ++d){      
  //     for (std::vector<int>::size_type n = 0; n < c0->neigh[F1][d].size(); ++n){
  //       int idn = c0->neigh[F1][d][n];
  //       Cell *cn = &Ctot[0][idn];
  //       double xn = cn->G.x[MV];
  //       double xLn = xn - cn->G.dx[MV]/2.;
  //       double xRn = xn + cn->G.dx[MV]/2.;

  //       Interface Int;
  //       Int.dim   = F1;
  //       Int.v     = 0;
  //       Int.x[F1] = (c0->G.x[F1] + cn->G.x[F1]) / 2.;
  //       Int.x[MV] = ( fmax(xL0,xLn) + fmin(xR0,xRn) )/2.;
  //       Int.dx[0] = fmax(0, fmin(xR0,xRn) - fmax(xL0,xLn));
  //       Int.computedA();
  //       dAtot += Int.dA;

  //       double xm = Int.x[MV]; 
  //       double dxmL  = xm - c0->G.cen[MV];
  //       double dxmR  = xm - cn->G.cen[MV];

  //       for (int q = 0; q < NUM_Q; ++q){
  //         double qL = c0->S.prim[q] + c0->grad_mv[q] * dxmL;
  //         double qR = cn->S.prim[q] + cn->grad_mv[q] * dxmR;
  //         double grad = (qR - qL)  / (cn->G.cen[F1] - c0->G.cen[F1]);
  //         grads[ind][q] = grad;
  //         grad_avg[q] += grad * Int.dA;
  //       }
  //       ind++;
  //     }
  //   }

  //   for (int q = 0; q < NUM_Q; ++q){
  //     grad_avg[q] /= dAtot;
  //     c0->grad_f1[q] = grad_avg[q]; // initialising for following minmod
  //   }

  //   // applying slope limiter
  //   for (int n = 0; n < n_neighs; ++n){
  //     for (int q = 0; q < NUM_Q; ++q){
  //       double tmp = c0->grad_f1[q];
  //       c0->grad_f1[q] = minmod(tmp, grads[n][q]);
  //     }
  //   }
  // }


#endif

void Grid::reconstructStates(int j, int i, int dim, int idn, Interface *Int){

  UNUSED(j);
  UNUSED(dim);
  UNUSED(idn);

  if (Int==NULL){ Int = &Itot[i]; }

  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_CONSTANT_
      Int->SL = Ctot[i  ].S;
      Int->SR = Ctot[i+1].S;
  #endif

  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_

    if (i == 0 or i == ntrack-2){
      Int->SL = Ctot[i  ].S;
      Int->SR = Ctot[i+1].S;
    } 
    else {
      Cell *cL  = &Ctot[i]; 
      Cell *cLL = &Ctot[max(i-1,0)]; 
      Cell *cR  = &Ctot[i+1]; 
      Cell *cRR = &Ctot[min(i+2,ntrack-1)];

      double gradL[NUM_Q], gradR[NUM_Q], gradN[NUM_Q];
      grad(*cLL, *cL, MV, gradL);
      grad(*cL,  *cR, MV, gradN);
      grad(*cR, *cRR, MV, gradR);
      for (int q = 0; q < NUM_Q; ++q){ 
        gradL[q] = minmod(gradL[q], gradN[q]);
        gradR[q] = minmod(gradR[q], gradN[q]);
        cL->grad_mv[q] = gradL[q];  // gradients stored in cells
        cR->grad_mv[q] = gradR[q];
      }

      Cell *c[] = {cL, cR};
      FluidState *IS[] = {&(Int->SL), &(Int->SR)};
      double *grad[] = {gradL, gradR};
      double xI = Int->x[dim];
      for (int s = 0; s < 2; ++s){
        for (int q = 0; q < NUM_Q; ++q){
          double xc = c[s]->G.cen[dim];
          double Sc = c[s]->S.prim[q];
          double g = grad[s][q];
          IS[s]->prim[q] = Sc + g * (xI - xc);
        }
        IS[s]->prim2cons(Int->x[x_]);
        IS[s]->state2flux(Int->x[x_]);
      }
    }


  //     Cell *cL = &Ctot[j][i]; 
  //     Cell *cR = &Ctot[0][idn];

  //     if (j == 0  
  //         or j == nde_nax[F1]-2
  //         or cL->neigh[F1][0].size() == 0
  //         or cR->neigh[F1][1].size() == 0){
  //       Int->SL = Ctot[j][i].S;
  //       Int->SR = Ctot[0][idn].S;
  //     }
  //     else {
  //       // recover gardients in moving direction and compute states at projected location
  //       // careful with the boundaries for LL and RR
  //       double xm = Int->x[MV];
  //       double xf = Int->x[F1];
  //       double dxmL  = xm - cL->G.cen[MV];
  //       double dxmR  = xm - cR->G.cen[MV];
  //       double dxfL  = xf - cL->G.cen[F1];
  //       double dxfR  = xf - cR->G.cen[F1];

  //       for (int q = 0; q < NUM_Q; ++q){
  //         double qL = cL->S.prim[q] + cL->grad_mv[q] * dxmL;
  //         double qR = cR->S.prim[q] + cR->grad_mv[q] * dxmR;
  //         double gradL = cL->grad_f1[q];
  //         double gradR = cR->grad_f1[q];
  //         Int->SL.prim[q] = qL + gradL * dxfL;
  //         Int->SR.prim[q] = qR + gradR * dxfR;
  //       }

  //       Int->SL.prim2cons(Int->x[x_]);
  //       Int->SR.prim2cons(Int->x[x_]);
  //       Int->SL.state2flux(Int->x[x_]);
  //       Int->SR.state2flux(Int->x[x_]);
  //     }
  //   }

  #endif

}


void Grid::computeNeighbors(bool print){

  UNUSED(print);

  for (int i = 1; i < ntrack-1; ++i){
    Cell *c = &Ctot[i];
    for (int d = 0; d < NUM_D; ++d){
      c->neigh[d][0].clear(); // resetting neighbors
      c->neigh[d][0].shrink_to_fit(); // resetting capacity
      c->neigh[d][1].clear();
      c->neigh[d][1].shrink_to_fit();

      c->neigh[d][0].push_back(Ctot[i-1].nde_id);
      c->neigh[d][1].push_back(Ctot[i+1].nde_id);
    }
  }

}


void Grid::regrid(){

  int j = 0;
  targetRegridVictims(j); // updating ismall and ibig
  for (int i = iLbnd+1; i <= iRbnd-1; ++i){  // only active cells
    int action = checkCellForRegrid(j, i);
    if (action != skip_){ 
      applyRegrid(j, i, action); 
    }
  }

}


void Grid::targetRegridVictims(int j){

  UNUSED(j);

  double minVal = 1.e20;
  double maxVal = 0.;
  for (int i = iLbnd+3; i <= iRbnd-3; ++i){   // not allowing edges to be victims
    Cell *c = &Ctot[i];

    double val = c->regridVal();
    if (val < 0) continue; 

    if (val < minVal){
      minVal = val;
      ismall = i;
    }
    if (val > maxVal){
      maxVal = val;
      ibig = i;
    }
  }

}


void Grid::applyRegrid(int j, int i, int action){


  // circular regrid: we grab the extra cell somewhere on the track (res stays the same)
  // And we avoid load-balancing issues
  if (action == split_ and ntrack < nde_nax[MV]-1){
    // merge smallest cell in track 
    int isplit = i;

    #if CIRC_REGRID_ == ENABLED_
      int is = ismall;
      if (i==is) exit(10);    // cell to split can't be the smallest cell in the fluid
      merge(j,is); // merging with smallest neighbour
      if (i > is) isplit--;
    #endif

    split(j,isplit);
  }

  if (action == merge_ and ntrack > 2*ngst+2){
    // we can just invert the idexes (i, is) -> (ib, i)
    merge(j,i); // merging with smallest neighbour

    #if CIRC_REGRID_ == ENABLED_
      int ib = ibig;
      if (i==ib) exit(10);    // cell to split can't be the smallest cell in the fluid
      if (ib > i) ib--;
      split(j,ib);      
    #endif

  }

  targetRegridVictims(j);  // updating ismall and ibig

}


void Grid::split(int j, int i){

  UNUSED(j);

  // freeing space on the right (so we can keep same i)
  std::copy_backward(&Ctot[i+1], &Ctot[ntrack],   &Ctot[ntrack+1]);
  std::copy_backward(&Itot[i],   &Itot[ntrack-1], &Itot[ntrack]  );
  ntrack++; // increasing the number of cells in track
  nact++; // increasing the number of cells in track
  iRbnd++;

  for (int ii = i; ii < ntrack; ++ii){
    int ind[]={ii};
    assignId(ind);
  }

  Cell *c  = &Ctot[i];
  Cell *cR = &Ctot[i+1];  // this should be an empty cell as we just freed it  


  // computing gradients
  // done before geometry update
  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_

    Cell *cL      = &Ctot[i-1]; 
    Cell *cR_old  = &Ctot[i+2];  // +2 because we freed up space on the right already

    // we're going to copy the conserved variables
    // c->S.prim2cons(c->G.cen[r_]); 
    // cL->S.prim2cons(cL->G.cen[r_]); 
    // cR_old->S.prim2cons(cR_old->G.cen[r_]); 

    double gradL[NUM_Q], gradR[NUM_Q], gradN[NUM_Q], gradC[NUM_Q];
    grad(*cL, *c, MV, gradL);
    grad(*c,  *cR_old, MV, gradR);
    grad(*cL,  *cR_old, MV, gradN);
    for (int q = 0; q < NUM_Q; ++q){ 
      gradL[q] = minmod(gradL[q], gradN[q]);
      gradR[q] = minmod(gradR[q], gradN[q]);
      gradC[q] = minmod(gradL[q], gradR[q]);
    }

  #endif

  // updating geometry
  double x_old = c->G.x[MV];
  double dl_old = c->G.dx[MV];
  // double cen_old = c->G.cen[MV];
  // double xL = x_old - dl_old/2.;
  // double xR = x_old + dl_old/2.;
  // double dlL = cen_old - xL;
  // double dlR = xR - cen_old;

  c->G.x[MV]  = x_old - dl_old/4.;
  cR->G.x[MV] = x_old + dl_old/4.;

  c->G.dx[MV]  = dl_old/2.;
  cR->G.dx[MV] = dl_old/2.;

  // c->G.x[MV]  = (xL+cen_old)/2.;
  // cR->G.x[MV] = (xR+cen_old)/2.;

  // c->G.dx[MV]  = dlL;
  // cR->G.dx[MV] = dlR;

  c->computeAllGeom();
  cR->computeAllGeom();

  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_CONSTANT_

    for (int q = 0; q < NUM_Q; ++q){
      cR->S.prim[q] = c->S.prim[q];
    }

  #endif
  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_

    double xL = c->G.cen[MV];
    double xR = cR->G.cen[MV];
    double xI = c->G.x[MV] + c->G.dx[MV]/2.;

    for (int q = 0; q < NUM_Q; ++q){
      double Sc = c->S.prim[q];
      double g = gradC[q];
      c->S.prim[q] = Sc + g * (xL - xI);
      cR->S.prim[q] = Sc + g * (xR - xI);
    }

  #endif

  // c->S.cons2prim(c->G.cen[r_]);
  // cR->S.cons2prim(cR->G.cen[r_]);

  c->S.prim2cons(c->G.cen[r_]);
  c->S.state2flux(c->G.cen[r_]);

  cR->S.prim2cons(cR->G.cen[r_]);
  cR->S.state2flux(cR->G.cen[r_]);

  Itot[i].x[MV]  = x_old;
  Itot[i].computedA();

}


void Grid::merge(int j, int i){

  UNUSED(j);

  Cell *c  = &Ctot[i];
  Cell *cL = &Ctot[i-1];
  Cell *cR = &Ctot[i+1];
  Cell *cVic;

  // identify neigbour victim
  int side;
  if ((cL->G.dx[MV] < cR->G.dx[MV] and i != ngst) or i == nact+ngst-1){
    side = left_;
    cVic = cL;
  } else {
    side = right_;
    cVic = cR;
  }

  c->S.prim2cons(c->G.cen[r_]);
  cVic->S.prim2cons(cVic->G.cen[r_]);
  // updating values
  double dV_loc = c->G.dV;
  double dV_vic = cVic->G.dV;
  for (int q = 0; q < NUM_Q; ++q){
    double Q_loc = c->S.cons[q];
    double Q_vic = cVic->S.cons[q];
    c->S.cons[q] = (Q_loc * dV_loc + Q_vic * dV_vic) / (dV_loc + dV_vic);
  }

  // updating geometry
  if (side == left_){
    c->G.x[MV] = ((cVic->G.x[MV] - cVic->G.dx[MV]/2.) + (c->G.x[MV] + c->G.dx[MV]/2.))/2.;
  }
  if (side == right_){
    c->G.x[MV] = ((cVic->G.x[MV] + cVic->G.dx[MV]/2.) + (c->G.x[MV] - c->G.dx[MV]/2.))/2.;
  }

  c->G.dx[MV] += cVic->G.dx[MV];
  c->computeAllGeom();

  c->S.cons2prim(c->G.cen[r_]);
  c->S.state2flux(c->G.cen[r_]);

  // shifting cells and interfaces around
  if (side==left_){
    std::copy(&Ctot[i], &Ctot[ntrack],   &Ctot[i-1]);
    std::copy(&Itot[i], &Itot[ntrack-1], &Itot[i-1]);
  }
  if (side==right_){
    std::copy(&Ctot[i+2], &Ctot[ntrack],   &Ctot[i+1]);
    std::copy(&Itot[i+1], &Itot[ntrack-1], &Itot[i]  );
  }
  ntrack--;  // not as many active cells
  nact--;
  iRbnd--;

  for (int ii = i-1; ii < ntrack; ++ii){
    int ind[]={ii};
    assignId(ind);
  }


}


void Grid::movDir_ComputeLambda(){

  int j=0;
  for (int i = 0; i < ntrack-1; ++i){
    reconstructStates(j,i,MV);
    Itot[i].computeLambda();
  }

}

void Grid::updateKinematics(int it, double t){

  for (int i = 0; i < ntrack-1; ++i){
    double v = VI * Itot[i].lS;
    double lfac = 1./sqrt(1.- v*v);
    Itot[i].v = v;
    Itot[i].lfac = lfac;
  }
  userKinematics(it, t);   // overriding with user preferences

}

void Grid::computeFluxes(){

  // flux in MV direction (looping over interfaces)
  for (int i = 0; i < ntrack-1; ++i){
    Itot[i].computeFlux();

    #if SHOCK_DETECTION_ == ENABLED_
      Cell *cL = &Ctot[i];
      Cell *cR = &Ctot[i+1];
      Itot[i].measureShock(cL, cR);
    #endif
      
  }
  for (int i = 1; i < ntrack-1; ++i){
    // can't do the edges
    
    Ctot[i].update_dt(MV, Itot[i-1], Itot[i]);
    for (int q = 0; q < NUM_Q; ++q){
      Ctot[i].flux[0][MV][q] = Itot[i-1].flux[q];
      Ctot[i].flux[1][MV][q] = Itot[i  ].flux[q];
    }
  }

}
  

double Grid::collect_dt(){

  double dt = 1.e15;
  double local_dt = 1.e15;
  for (int i = 0; i < ntrack; ++i){
    local_dt = fmin(local_dt, Ctot[i].dt_loc);
    Ctot[i].dt_loc = 1.e15;  // resetting dt_loc for next update
  }

  // // taking minimum dt accross all processes
  // MPI_Barrier(MPI_COMM_WORLD); // wait until all processes have caught up
  // MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  dt = local_dt;
  return dt;

}


void Grid::update(double dt){

  for (int i = 0; i < ntrack-1; ++i){
    Itot[i].move(dt);
  }
  // do not update border cells because can lead to non-physical states
  for (int i = iLbnd+1; i <= iRbnd-1; ++i){
    double xL = Itot[i-1].x[MV];
    double xR = Itot[i].x[MV];
    Ctot[i].update(dt,xL,xR);
  }

}


void Grid::copyState0(){
  // Copies the current state of the grid into S0 for higher order time-stepping

  for (int i = 0; i < ntrack; ++i){
    Cell *c = &Ctot[i];
    c->S0 = c->S;
    c->G0 = c->G;
  }

  for (int i = 0; i < ntrack-1; ++i){
    Interface *I = &Itot[i];
    for (int d = 0; d < NUM_D; ++d) I->x0[d] = I->x[d];
    I->v0 = I->v;
  }

}


void Grid::CellGeomFromInterfacePos(){
  // pos in F1 needs to already be updated

  for (int i = 1; i < ntrack-1; ++i){
    Cell *c = &Ctot[i];
    double xL = Itot[i-1].x[MV];
    double xR = Itot[i].x[MV];
    c->move(xL, xR);
  }

}


void Grid::interfaceGeomFromCellPos(){

  for (int i = 0; i < ntrack-1; ++i){
    if (i != iLbnd){
      Itot[i].x[MV] = Ctot[i].G.x[MV] + Ctot[i].G.dx[MV]/2.;
    }
    else{
      Itot[i].x[MV] = Ctot[i+1].G.x[MV] - Ctot[i+1].G.dx[MV]/2.;
    }
    Itot[i].computedA();
  }

}


void Grid::interfaceGeomFromCellPos(int j){

  UNUSED(j);
  for (int i = 0; i < ntrack-1; ++i){
    Itot[i].x[MV] = Ctot[i+1].G.x[MV] - Ctot[i+1].G.dx[MV]/2.;
    Itot[i].computedA();
  }

}


void Grid::destruct(){

  delete_array_1d(Ctot);
  delete_array_1d(Itot);

}


void Grid::apply(void (Cell::*func)()){

  // usage: grid.apply(&Cell::func_name);
  
  for (int i = 0; i < nde_nax[MV]; ++i){
    (Ctot[i].*func)();
  }

}


// overloading apply for FluidState
// void Grid::apply(void (FluidState::*func)(), bool noborder){

//   int jL, jR;
//   if (noborder){ jL = 1; jR = nde_nax[F1]-1; }
//   else         { jL = 0; jR = nde_nax[F1];   }

//   for (int j = jL; j < jR; ++j){

//     int iL, iR;
//     if (noborder){ iL = 1; iR = ntrack[j]-1; }
//     else         { iL = 0; iR = ntrack[j];   }

//     for (int i = iL; i < iR; ++i){
//       (Ctot[j][i].S.*func)();
//     }
//   }
    
// }


void Grid::prim2cons(){

  int iL = 0;
  int iR = ntrack;
  for (int i = iL; i < iR; ++i){
    double r = Ctot[i].G.cen[r_];
    Ctot[i].S.prim2cons(r);
  }

}


void Grid::state2flux(){

  int iL = 0;
  int iR = ntrack;
  for (int i = iL; i < iR; ++i){
    double r = Ctot[i].G.cen[r_];
    Ctot[i].S.state2flux(r);
  }

}


void Grid::v2u(){

  // converting velocities to 4-velocities:
  for (int i = 0; i < ntrack; ++i){
    double v   = 0;
    double vv[NUM_D];
    for (int d = 0; d < NUM_D; ++d){
      vv[d] = Ctot[i].S.prim[VV1+d];
      v += vv[d]*vv[d];
    }
    v = sqrt(v);
    double lfac = 1./sqrt(1.-v*v);
    for (int d = 0; d < NUM_D; ++d){
      Ctot[i].S.prim[UU1+d] *= lfac;
    }
  }

}




void mpi_distribute(Grid *grid){

  int ngst = grid->ngst;
  int size  = grid->nde_nax[MV]-2*ngst;
  std::copy_n(&(grid->Cinit[0]), size, &(grid->C[0]));
  delete_array_1d(grid->Cinit);

  // storing index information in cell
  for (int i = 0; i < grid->nde_nax[MV]; ++i){
    grid->Ctot[i].nde_id = i;
    grid->Ctot[i].nde_ind[0] = i;
  }

}







