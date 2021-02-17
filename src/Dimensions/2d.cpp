/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-02-17 11:40:38
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


template <> Interface** array_2d_nogst<Interface>(Interface** in, int nj, const int ngst){

  Interface** pp;
  int nj_new = nj-2*ngst;
  try { pp = new Interface*[nj_new]; }
  catch (bad_alloc&) {
    cout<<"arr_2d: memory allocation error, pp = new T*["<<nj_new<<"]"<<endl;
    exit(1);
  }
  for (int j=0; j<nj_new; ++j) { pp[j] = in[j+ngst] + (ngst-1); }
  return pp;

}

static int splitGrid(int ncell, int ngst, int *origin){
  /*
  int ncell:   total number of cells (non including ghost boundaries)
  int ngst:    size of ghost regions.
  int *origin:  coordinates of C[0][0] in the simulation domain
  */
  int base = ncell/worldsize;
  int rest = ncell%worldsize;
  int nde_nax = base;
  if (worldrank<rest) nde_nax++;
  nde_nax+=2*ngst;

  origin[MV] = 0;
  origin[F1] = worldrank*base;
  if (worldrank <rest) origin[F1]+=worldrank;
  else origin[F1]+=rest;

  return nde_nax;

}

void Grid::initialise(s_par par){

  for (int d = 0; d < NUM_D; ++d){
    ncell[d] = par.ncell[d];
    ngst     = par.ngst;
    nax[d]   = ncell[d]+2*ngst;
  }
  nax[MV]    = par.nmax+2*ngst;  // overriding with max number of cells
  nsimu = 1; for (int d = 0; d < NUM_D; ++d){ nsimu *= ncell[d]; }
  nde_nax[MV]   = nax[MV];
  nde_nax[F1]   = splitGrid(ncell[F1], ngst, origin);
  nde_ncell[MV] = ncell[MV];  // max number of active cells in mov dim
  nde_ncell[F1] = nde_nax[F1]-2*ngst;
  nde_ntot      = nde_nax[F1] * nde_nax[MV];

  // active cells and ghost cells
  nact  = new int[nde_nax[F1]];
  for (int j = 0; j < nde_nax[F1]; ++j){nact[j] = nde_ncell[MV];}
  ntrack  = new int[nde_nax[F1]];
  for (int j = 0; j < nde_nax[F1]; ++j){ntrack[j] = nact[j]+2*ngst;}
  jLbnd = ngst-1;
  jRbnd = nde_nax[F1]-ngst;
  iLbnd = new int[nde_nax[F1]];
  iRbnd = new int[nde_nax[F1]];
  for (int j = 0; j < nde_nax[F1]; ++j){
    iLbnd[j] = ngst-1;
    iRbnd[j] = nde_ncell[MV]+ngst;
  }
  ibig   = new int[nde_nax[F1]];
  ismall = new int[nde_nax[F1]];
  
  Cinit = array_2d<Cell>(ncell[F1],nax[MV]-2*ngst); // Cinit need to have the empty cells...
  Ctot  = array_2d<Cell>(nde_nax[F1],nde_nax[MV]);                       // ...for restart
  Itot  = array_2d<Interface>(nde_nax[F1],nde_nax[MV]-1);
  I     = array_2d_nogst<Interface>(Itot, nde_nax[F1], ngst); // removes ghost to ghost
  C     = array_2d_nogst<Cell>(Ctot, nde_nax[F1], ngst);     

}


void Grid::assignId(int ind[NUM_D]){

  Cell *c = &Ctot[ind[0]][ind[1]];
  c->nde_id = ind[0]*nde_nax[MV] + ind[1];
  for (int d = 0; d < NUM_D; ++d){ c->nde_ind[d] = ind[d]; }

}


void Grid::mpi_exchangeGhostTracks(){

  int rL = worldrank-1;
  int rR = worldrank+1;
  MPI_Datatype cell_mpi = {0}; 
  generate_mpi_cell( &cell_mpi );
  MPI_Barrier(MPI_COMM_WORLD);

  // for interfaces I just need to exchange the positions. That can be done in the cells
  // directly or I can just exchange a vector of positions

  // sending to lower node, receiving from higher node
  for (int j = 0; j < ngst; ++j){
    int jout = j+ngst;
    int jin  = jRbnd+j;
    int nout = nact[jout]+2*ngst; // number of cells to receive
    int nin  = 1;                 // number of cells to send (default is 1)
    int tag1 = j;
    int tag2 = ngst+j;

    if (worldrank==0){
      MPI_Recv(&nin , 1, MPI_INT, rR, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
    } else if (worldrank==worldsize-1){
      MPI_Send(&nout, 1, MPI_INT, rL, tag1, MPI_COMM_WORLD);      
    } else {
      MPI_Sendrecv(&nout, 1, MPI_INT, rL, tag1, 
                   &nin , 1, MPI_INT, rR, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
    }

    // updating ntrack, nact, iRbnd (iLbnd doesn't need to be updated)
    if (worldrank!=worldsize-1){
      ntrack[jin] = nin;
      nact[jin]  = nin-2*ngst;
      iRbnd[jin] = nin-ngst; 
    }

    s_cell *SCout = new s_cell[nout];
    s_cell *SCin  = new s_cell[nin];
    if (worldrank==0){
      MPI_Recv(SCin , nin , cell_mpi, rR, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }    
    } 
    else if (worldrank==worldsize-1){
      for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
      MPI_Send(SCout, nout, cell_mpi, rL, tag2, MPI_COMM_WORLD);      
    } 
    else {
      for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
      MPI_Sendrecv(SCout, nout, cell_mpi, rL, tag2, 
                   SCin , nin , cell_mpi, rR, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }          
    }
    for (int i = 0; i < nin; ++i){ delete &SCin[i]; }
    for (int i = 0; i < nout; ++i){ delete &SCout[i]; }
    delete[](SCin); delete[](SCout);
  }

  // sending to higher node, receiving from lower node
  for (int j = 0; j < ngst; ++j){
    int jout = jRbnd-ngst+j;
    int jin  = j;
    int nout = nact[jout]+2*ngst; // number of cells to receive
    int nin  = 1;                 // number of cells to send (default is 1)
    int tag1 = j;
    int tag2 = ngst+j;

    if (worldrank==worldsize-1){
      MPI_Recv(&nin , 1, MPI_INT, rL, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
    } else if (worldrank==0){
      MPI_Send(&nout, 1, MPI_INT, rR, tag1, MPI_COMM_WORLD);      
    } else {
      MPI_Sendrecv(&nout, 1, MPI_INT, rR, tag1, 
                   &nin , 1, MPI_INT, rL, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // updating ntrack, nact, iRbnd (iLbnd doesn't need to be updated)
    if (worldrank!=0){
      ntrack[jin] = nin;
      nact[jin]  = nin-2*ngst;
      iRbnd[jin] = nin-ngst;
    }

    s_cell *SCout = new s_cell[nout];
    s_cell *SCin  = new s_cell[nin];
    if (worldrank==worldsize-1){
      MPI_Recv(SCin , nin , cell_mpi, rL, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }
    } 
    else if (worldrank==0){
      for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
      MPI_Send(SCout, nout, cell_mpi, rR, tag2, MPI_COMM_WORLD);      
    } 
    else {
      for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
      MPI_Sendrecv(SCout, nout, cell_mpi, rR, tag2, 
                   SCin , nin , cell_mpi, rL, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }          
    }
    for (int i = 0; i < nin; ++i){ delete &SCin[i]; }
    for (int i = 0; i < nout; ++i){ delete &SCout[i]; }
    delete[](SCin); delete[](SCout);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int j = 0; j < ngst; ++j){ // looping over ghost tracks
    for (int i = 0; i < nde_nax[MV]; ++i){

      // updating interface positions from communicated cells
      if (i != nde_nax[MV]-1){ 
        int jn = nde_nax[F1]-1-j;
        Itot[j][i].x[MV]  = Ctot[j][i].G.x[MV] + Ctot[j][i].G.dx[MV]/2.;
        Itot[jn][i].x[MV] = Ctot[jn][i].G.x[MV] + Ctot[jn][i].G.dx[MV]/2.;
        Itot[j][i].x[F1]  = Ctot[j][i].G.x[F1];
        Itot[jn][i].x[F1] = Ctot[jn][i].G.x[F1];
      }

      // updating IDs
      int indL[] = {j,i};
      int indR[] = {nde_nax[F1]-1-j,i};
      assignId(indL);
      assignId(indR);
    }
  }

}


void Grid::updateGhosts(int it, double t){

  if (worldsize>1) { mpi_exchangeGhostTracks(); }

  // outer boundaries (OUTFLOW)
  if (worldrank==0){
    for (int j = 0; j <= jLbnd; ++j){
      ntrack[j] = ntrack[jLbnd+1];
      nact[j] = nact[jLbnd+1];
      iRbnd[j] = iRbnd[jLbnd+1];
      iLbnd[j] = iLbnd[jLbnd+1];
      std::copy_n(&Ctot[jLbnd+1][0], ntrack[jLbnd+1],   &Ctot[j][0]);
      std::copy_n(&Itot[jLbnd+1][0], ntrack[jLbnd+1]-1, &Itot[j][0]);

      // updating ghost positions, ids and indexes
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1]; 
          Itot[j][i].computedA();
        }
      }
    }
  }
  if (worldrank==worldsize-1){
    for (int j = jRbnd; j < nde_nax[F1]; ++j){
      ntrack[j] = ntrack[jRbnd-1];
      nact[j] = nact[jRbnd-1];
      iRbnd[j] = iRbnd[jRbnd-1];
      iLbnd[j] = iLbnd[jRbnd-1];
      std::copy_n(&Ctot[jRbnd-1][0], ntrack[jRbnd-1]  , &Ctot[j][0]);
      std::copy_n(&Itot[jRbnd-1][0], ntrack[jRbnd-1]-1, &Itot[j][0]);

      // udating ghost positions
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
          Itot[j][i].computedA();
        }
      }
    }
  }
  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j)
  {
    for (int i = 0; i <= iLbnd[j]; ++i){
      Ctot[j][i] = Ctot[j][iLbnd[j]+1];
      Ctot[j][i].G.x[MV] -= (iLbnd[j]-i+1) * Ctot[j][iLbnd[j]+1].G.dx[MV];
      Ctot[j][i].computeAllGeom();
      int ind[] = {j,i};
      assignId(ind);
    }
    for (int i = 0; i <= iLbnd[j]-1; ++i){
      Itot[j][i] = Itot[j][iLbnd[j]+1];
      Itot[j][i].x[MV] -= (iLbnd[j]-i+1) * Ctot[j][iLbnd[j]+1].G.dx[MV];
      Itot[j][i].computedA();
    }

    for (int i = iRbnd[j]; i < ntrack[j]; ++i){
      Ctot[j][i] = Ctot[j][iRbnd[j]-1];
      Ctot[j][i].G.x[MV] += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dx[MV];
      Ctot[j][i].computeAllGeom();
      int ind[] = {j,i};
      assignId(ind);
    }
    for (int i = iRbnd[j]; i < ntrack[j]-1; ++i){
      Itot[j][i] = Itot[j][iRbnd[j]-1];
      Itot[j][i].x[MV] += (i-iRbnd[j]) * Ctot[j][iRbnd[j]-1].G.dx[MV];
      Itot[j][i].computedA();
    }
  }
  userBoundaries(it, t); // overiding with user-specific boundary conditions

}

#if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_

  static double minmod(double a, double b){

    if (fabs(a) < fabs(b) && a*b > 0){ return(a); }
    if (fabs(b) < fabs(a) && a*b > 0){ return(b); }
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

  void Grid::computeTransGradients(int j, int i){

    Cell *c0 = &Ctot[j][i];
    double x0 = c0->G.x[MV];
    double xL0 = x0 - c0->G.dx[MV]/2.;
    double xR0 = x0 + c0->G.dx[MV]/2.;
    const int n_neighs = (int) (c0->neigh[F1][0].size() + c0->neigh[F1][1].size());
    double grads[n_neighs][NUM_Q];
    double grad_avg[NUM_Q] = {0};
    double dAtot = 0;

    int ind = 0;
    for (int d = 0; d < 2; ++d){      
      for (std::vector<int>::size_type n = 0; n < c0->neigh[F1][d].size(); ++n){
        int idn = c0->neigh[F1][d][n];
        Cell *cn = &Ctot[0][idn];
        double xn = cn->G.x[MV];
        double xLn = xn - cn->G.dx[MV]/2.;
        double xRn = xn + cn->G.dx[MV]/2.;

        double xf1I = c0->G.x[F1];
        if (d == 0){
          xf1I -= c0->G.dx[F1]/2.;
        } else if (d==0){
          xf1I += c0->G.dx[F1]/2.;
        }
        Interface Int;
        Int.dim   = F1;
        Int.v     = 0;
        // Int.x[F1] = (c0->G.x[F1] + cn->G.x[F1]) / 2.;
        Int.x[F1] = xf1I;
        Int.x[MV] = ( fmax(xL0,xLn) + fmin(xR0,xRn) )/2.;
        Int.dx[0] = fmax(0, fmin(xR0,xRn) - fmax(xL0,xLn));
        Int.computedA();
        dAtot += Int.dA;

        double xm = Int.x[MV]; 
        double dxmL  = xm - c0->G.cen[MV];
        double dxmR  = xm - cn->G.cen[MV];

        for (int q = 0; q < NUM_Q; ++q){
          double qL = c0->S.prim[q] + c0->grad_mv[q] * dxmL;
          double qR = cn->S.prim[q] + cn->grad_mv[q] * dxmR;
          double grad = (qR - qL)  / (cn->G.cen[F1] - c0->G.cen[F1]);
          grads[ind][q] = grad;
          grad_avg[q] += grad * Int.dA;
        }
        ind++;
      }
    }

    for (int q = 0; q < NUM_Q; ++q){
      grad_avg[q] /= dAtot;
      c0->grad_f1[q] = grad_avg[q]; // initialising for following minmod
    }

    // applying slope limiter
    for (int n = 0; n < n_neighs; ++n){
      for (int q = 0; q < NUM_Q; ++q){
        double tmp = c0->grad_f1[q];
        c0->grad_f1[q] = minmod(tmp, grads[n][q]);
      }
    }
  }


#endif

void Grid::reconstructStates(int j, int i, int dim, int idn, Interface *Int){

  if (Int==NULL){ Int = &Itot[j][i]; }

  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_CONSTANT_
    if (dim == MV){
      Int->SL = Ctot[j][i  ].S;
      Int->SR = Ctot[j][i+1].S;
      UNUSED(idn);
    } else {
      Int->SL = Ctot[j][i].S;
      Int->SR = Ctot[0][idn].S;
    }
  #endif

  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_
    if (dim == MV){

      if (i == 0 or i == ntrack[j]-2){
        Int->SL = Ctot[j][i  ].S;
        Int->SR = Ctot[j][i+1].S;
      } 
      else {
        Cell *cL  = &Ctot[j][i]; 
        Cell *cLL = &Ctot[j][max(i-1,0)]; 
        Cell *cR  = &Ctot[j][i+1]; 
        Cell *cRR = &Ctot[j][min(i+2,ntrack[j]-1)];

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
      UNUSED(idn);

    } else {

      Cell *cL = &Ctot[j][i]; 
      Cell *cR = &Ctot[0][idn];

      if (j == 0  
          or j == nde_nax[F1]-2
          or cL->neigh[F1][0].size() == 0
          or cR->neigh[F1][1].size() == 0){
        Int->SL = Ctot[j][i].S;
        Int->SR = Ctot[0][idn].S;
      }
      else {
        // recover gardients in moving direction and compute states at projected location
        // careful with the boundaries for LL and RR
        double xm = Int->x[MV];
        double xf = Int->x[F1];
        double dxmL  = xm - cL->G.cen[MV];
        double dxmR  = xm - cR->G.cen[MV];
        double dxfL  = xf - cL->G.cen[F1];
        double dxfR  = xf - cR->G.cen[F1];

        for (int q = 0; q < NUM_Q; ++q){
          double qL = cL->S.prim[q] + cL->grad_mv[q] * dxmL;
          double qR = cR->S.prim[q] + cR->grad_mv[q] * dxmR;
          double gradL = cL->grad_f1[q];
          double gradR = cR->grad_f1[q];
          Int->SL.prim[q] = qL + gradL * dxfL;
          Int->SR.prim[q] = qR + gradR * dxfR;
        }

        Int->SL.prim2cons(Int->x[x_]);
        Int->SR.prim2cons(Int->x[x_]);
        Int->SL.state2flux(Int->x[x_]);
        Int->SR.state2flux(Int->x[x_]);
      }
    }

  #endif

}


void Grid::computeNeighbors(bool print){

  for (int j = 0; j < nde_nax[F1]; ++j){ // can't do it for edges
    int im = 0;
    int ip = 0;
    for (int i = 1; i < ntrack[j]-1; ++i){
      Cell *c = &Ctot[j][i];
      for (int d = 0; d < NUM_D; ++d){
        c->neigh[d][0].clear(); // resetting neighbors
        c->neigh[d][0].shrink_to_fit(); // resetting capacity
        c->neigh[d][1].clear();
        c->neigh[d][1].shrink_to_fit();

        if (d == MV){
          c->neigh[d][0].push_back(Ctot[j][i-1].nde_id);
          c->neigh[d][1].push_back(Ctot[j][i+1].nde_id);
        } else {
          double xjL = c->G.x[MV] - c->G.dx[MV]/2.;
          double xjR = c->G.x[MV] + c->G.dx[MV]/2.;

          if (j!=0){
            double xm = Itot[j-1][im].x[MV];
            if (xm > xjL and im > 0){ 
              im--; xm = Itot[j-1][im].x[MV];
            }
            while (xm < xjR){
              c->neigh[d][0].push_back(Ctot[j-1][im+1].nde_id);
              if (im >= ntrack[j-1]-2) break;
              im++; xm = Itot[j-1][im].x[MV];
            }            
          }

          if (j!=nde_nax[F1]-1){
            double xp = Itot[j+1][ip].x[MV];
            if (xp > xjL and ip > 0){ 
              ip--; xp = Itot[j+1][ip].x[MV];
              } 
            while (xp < xjR){
              c->neigh[d][1].push_back(Ctot[j+1][ip+1].nde_id);
              if (ip >= ntrack[j+1]-2) break;
              ip++; xp = Itot[j+1][ip].x[MV];
            }
          }
        }
      }
    }
  }

  if (print){
    printf("%d %d\n", MV, F1);
    for (int j = 1; j < nde_nax[F1]-1; ++j){
      for (int i = 1; i < ntrack[j]-1; ++i){
        for (int d = 0; d < NUM_D; ++d){
          for (int s = 0; s < 2; ++s)
          {
            Cell *c  = &Ctot[j][i];
            printf("%d %d, %d %d, %le: ", j, i, d, s, c->G.x[MV]);
            for (std::vector<int>::size_type n = 0; n < c->neigh[d][s].size(); ++n){
              int id  = c->neigh[d][s][n];
              Cell *cn = &Ctot[0][id];
              int jn  = cn->nde_ind[0];
              int in  = cn->nde_ind[1];
              printf(" %d %d %le %le", jn, in, cn->G.x[MV], cn->G.dx[MV]);
            }
            printf("\n");  
          }
        }
      }
    } 
    exit(9);
  }

}


void Grid::regrid(){

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    targetRegridVictims(j); // updating ismall and ibig
    for (int i = iLbnd[j]+1; i <= iRbnd[j]-1; ++i){  // only active cells
      int action = checkCellForRegrid(j, i);
      if (action != skip_){ 
        applyRegrid(j, i, action); 
      }
    }
  }

}


void Grid::targetRegridVictims(int j){

  double minVal = 1.e20;
  double maxVal = 0.;
  for (int i = iLbnd[j]+3; i <= iRbnd[j]-3; ++i){   // not allowing edges to be victims
    Cell *c = &Ctot[j][i];

    double val = c->regridVal();
    if (val < 0) continue; 

    if (val < minVal){
      minVal = val;
      ismall[j] = i;
    }
    if (val > maxVal){
      maxVal = val;
      ibig[j] = i;
    }
  }

}


void Grid::applyRegrid(int j, int i, int action){

  // printf("%d\n", action);

  // circular regrid: we grab the extra cell somewhere on the track (res stays the same)
  // And we avoid load-balancing issues
  if (action == split_ and ntrack[j] < nde_nax[MV]-1){
    // merge smallest cell in track 
    int isplit = i;

    #if CIRC_REGRID_ == ENABLED_
      int is = ismall[j];
      if (i==is) exit(10);    // cell to split can't be the smallest cell in the fluid
      merge(j,is); // merging with smallest neighbour
      if (i > is) isplit--;
    #endif

    split(j,isplit);
  }

  if (action == merge_ and ntrack[j] > 2*ngst+2){
    // we can just invert the idexes (i, is) -> (ib, i)
    merge(j,i); // merging with smallest neighbour

    #if CIRC_REGRID_ == ENABLED_
      int ib = ibig[j];
      if (i==ib) exit(10);    // cell to split can't be the smallest cell in the fluid
      if (ib > i) ib--;
      split(j,ib);      
    #endif

  }

  targetRegridVictims(j);  // updating ismall and ibig

}


void Grid::split(int j, int i){


  // freeing space on the right (so we can keep same i)
  std::copy_backward(&Ctot[j][i+1], &Ctot[j][ntrack[j]],   &Ctot[j][ntrack[j]+1]);
  std::copy_backward(&Itot[j][i],   &Itot[j][ntrack[j]-1], &Itot[j][ntrack[j]]  );
  ntrack[j]++; // increasing the number of cells in track
  nact[j]++; // increasing the number of cells in track
  iRbnd[j]++;

  for (int ii = i; ii < ntrack[j]; ++ii){
    int ind[]={j,ii};
    assignId(ind);
  }

  Cell *c  = &Ctot[j][i];
  Cell *cR = &Ctot[j][i+1];  // this should be an empty cell as we just freed it  

  // computing gradients
  // done before geometry update
  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_LINEAR_

    Cell *cL      = &Ctot[j][i-1]; 
    Cell *cR_old  = &Ctot[j][i+2];  // +2 because we freed up space on the right already

    double gradL[NUM_Q], gradR[NUM_Q], gradC[NUM_Q];
    grad(*cL, *c, MV, gradL);
    grad(*c,  *cR_old, MV, gradR);
    for (int q = 0; q < NUM_Q; ++q){ 
      gradC[q] = minmod(gradL[q], gradR[q]);
    }

  #endif

  // updating geometry
  double x_old = c->G.x[MV];
  double dl_old = c->G.dx[MV];

  c->G.x[MV]  = x_old - dl_old/4.;
  cR->G.x[MV] = x_old + dl_old/4.;
  cR->G.x[F1] = c->G.x[F1];

  c->G.dx[MV]  = dl_old/2.;
  cR->G.dx[MV] = dl_old/2.;
  cR->G.dx[F1] = c->G.dx[F1];

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

  c->S.prim2cons(c->G.x[r_]);
  c->S.state2flux(c->G.x[r_]);

  cR->S.prim2cons(cR->G.x[r_]);
  cR->S.state2flux(cR->G.x[r_]);

  Itot[j][i].x[MV]  = x_old;
  Itot[j][i].x[F1]  = c->G.x[F1];
  Itot[j][i].dx[0] = c->G.dx[F1];
  Itot[j][i].computedA();

}


void Grid::merge(int j, int i){

  Cell *c  = &Ctot[j][i];
  Cell *cL = &Ctot[j][i-1];
  Cell *cR = &Ctot[j][i+1];
  Cell *cVic;

  // identify neigbour victim
  int side;
  if ((cL->G.dx[MV] < cR->G.dx[MV] and i != ngst) or i == nact[j]+ngst-1){
    side = left_;
    cVic = cL;
  } else {
    side = right_;
    cVic = cR;
  }

  c->S.prim2cons(c->G.x[r_]);
  cVic->S.prim2cons(cVic->G.x[r_]);
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

  c->S.cons2prim(c->G.x[r_]);
  c->S.state2flux(c->G.x[r_]);

  // shifting cells and interfaces around
  if (side==left_){
    std::copy(&Ctot[j][i], &Ctot[j][ntrack[j]],   &Ctot[j][i-1]);
    std::copy(&Itot[j][i], &Itot[j][ntrack[j]-1], &Itot[j][i-1]);
  }
  if (side==right_){
    std::copy(&Ctot[j][i+2], &Ctot[j][ntrack[j]],   &Ctot[j][i+1]);
    std::copy(&Itot[j][i+1], &Itot[j][ntrack[j]-1], &Itot[j][i]  );
  }
  ntrack[j]--;  // not as many active cells
  nact[j]--;
  iRbnd[j]--;

  for (int ii = i-1; ii < ntrack[j]; ++ii){
    int ind[]={j,ii};
    assignId(ind);
  }


}


void Grid::movDir_ComputeLambda(){

  #pragma omp parallel for default(shared)
  // this needs to be done on all j tracks because we use the gradients later on
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      reconstructStates(j,i,MV);
      Itot[j][i].computeLambda();
    }
  }

}

void Grid::updateKinematics(int it, double t){

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      double v = VI * Itot[j][i].lS;
      double lfac = 1./sqrt(1.- v*v);
      Itot[j][i].v = v;
      Itot[j][i].lfac = lfac;
    }
  }  
  userKinematics(it, t);   // overriding with user preferences

}

void Grid::computeFluxes(){

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){

    // flux in MV direction (looping over interfaces)
    for (int i = 0; i < ntrack[j]-1; ++i){
      Itot[j][i].computeFlux();

      #if SHOCK_DETECTION_ == ENABLED_
        Cell *cL = &Ctot[j][i];
        Cell *cR = &Ctot[j][i+1];
        Itot[j][i].measureShock(cL, cR);
      #endif
        
    }
    for (int i = 1; i < ntrack[j]-1; ++i){
      // can't do the edges
      
      Ctot[j][i].update_dt(MV, Itot[j][i-1], Itot[j][i]);
      for (int q = 0; q < NUM_Q; ++q){
        Ctot[j][i].flux[0][MV][q] = Itot[j][i-1].flux[q];
        Ctot[j][i].flux[1][MV][q] = Itot[j][i  ].flux[q];
      }
    }
  }
  // flux in F1 direction (building the interfaces from neighbor ids)
  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]-1; ++j){
    // double jpos = ( Ctot[j][iLbnd[j]+1].G.x[F1] + Ctot[j+1][iLbnd[j+1]+1].G.x[F1] )/2.;
    double jpos = Ctot[j][iLbnd[j]+1].G.x[F1] + Ctot[j][iLbnd[j]+1].G.dx[F1]/2.;

    // resetting fluxes
    for (int i = 0; i < ntrack[j]; ++i){
      for (int q = 0; q < NUM_Q; ++q){ Ctot[j][i].flux[1][F1][q] = 0.; }
    }
    for (int i = 0; i < ntrack[j+1]; ++i){
      for (int q = 0; q < NUM_Q; ++q){ Ctot[j+1][i].flux[0][F1][q] = 0.; }
    }
    for (int i = 0; i < ntrack[j]; ++i){

      Cell *c0 = &Ctot[j][i];
      double x0 = c0->G.x[MV];
      double xL0 = x0 - c0->G.dx[MV]/2.;
      double xR0 = x0 + c0->G.dx[MV]/2.;

      computeTransGradients(j,i);

      // We reconstruct on the (+) side of the track (or R-side)
      // this will be taken into account in the choice of gradients of reconstructStates()
      for (std::vector<int>::size_type n = 0; n < c0->neigh[F1][1].size(); ++n){
        int idn = c0->neigh[F1][1][n];
        Cell *cn = &Ctot[0][idn];
        double xn = cn->G.x[MV];
        double xLn = xn - cn->G.dx[MV]/2.;
        double xRn = xn + cn->G.dx[MV]/2.;

        Interface Int;
        Int.dim   = F1;
        Int.v     = 0;
        Int.x[F1] = jpos;
        Int.x[MV] = ( fmax(xL0,xLn) + fmin(xR0,xRn) )/2.;
        Int.dx[0] = fmax(0, fmin(xR0,xRn) - fmax(xL0,xLn));
        Int.computedA();

        reconstructStates(j,i,F1,idn,&Int);
        Int.computeLambda();
        Int.computeFlux();

        #if SHOCK_DETECTION_ == ENABLED_
          // printf("%le %le %d %d\n", c0->S.prim[PPP], cn->S.prim[PPP], c0->nde_ind[1], cn->nde_ind[1]);
          Int.measureShock(c0, cn);
        #endif

        c0->update_dt(F1, Int.lL);
        cn->update_dt(F1, Int.lR);
        for (int q = 0; q < NUM_Q; ++q){
          c0->flux[1][F1][q] += Int.flux[q];
          cn->flux[0][F1][q] += Int.flux[q];
        }
      }
    }  
  }
}
  

double Grid::collect_dt(){

  double dt = 1.e15;
  double local_dt = 1.e15;
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]; ++i){
      local_dt = fmin(local_dt, Ctot[j][i].dt_loc);
      Ctot[j][i].dt_loc = 1.e15;  // resetting dt_loc for next update
    }
  }

  // taking minimum dt accross all processes
  MPI_Barrier(MPI_COMM_WORLD); // wait until all processes have caught up
  MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return dt;

}


void Grid::update(double dt){

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      Itot[j][i].move(dt);
    }
  }
  // do not update border cells because can lead to non-physical states
  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]+1; i <= iRbnd[j]-1; ++i){
      double xL = Itot[j][i-1].x[MV];
      double xR = Itot[j][i].x[MV];
      double vL = Itot[j][i-1].v;
      double vR = Itot[j][i].v;
      if (xR < xL) printf("vel %d %d %le %le %le %le\n", j, i, vL, vR, xL, xR);
      // printf("%d %d\n", j, i);
      Ctot[j][i].update(dt,xL,xR);
    }
  }

}


void Grid::copyState0(){
  // Copies the current state of the grid into S0 for higher order time-stepping

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]; ++i){
      Cell *c = &Ctot[j][i];
      c->S0 = c->S;
      c->G0 = c->G;
    }
  }  

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      Interface *I = &Itot[j][i];
      for (int d = 0; d < NUM_D; ++d) I->x0[d] = I->x[d];
      I->v0 = I->v;
    }
  }  

}


void Grid::CellGeomFromInterfacePos(){
  // pos in F1 needs to already be updated

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 1; i < ntrack[j]-1; ++i){
      Cell *c = &Ctot[j][i];
      double xL = Itot[j][i-1].x[MV];
      double xR = Itot[j][i].x[MV];
      c->move(xL, xR);
    }
  }

}


void Grid::interfaceGeomFromCellPos(){

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    double xj = Ctot[j][iLbnd[j]+1].G.x[F1];
    for (int i = 0; i < ntrack[j]-1; ++i){
      if (i != iLbnd[j]){
        Itot[j][i].x[MV] = Ctot[j][i].G.x[MV] + Ctot[j][i].G.dx[MV]/2.;
        Itot[j][i].dx[0] = Ctot[j][i].G.dx[F1];
      }
      else{
        Itot[j][i].x[MV] = Ctot[j][i+1].G.x[MV] - Ctot[j][i+1].G.dx[MV]/2.;
        Itot[j][i].dx[0] = Ctot[j][i+1].G.dx[F1];
      }
      Itot[j][i].x[F1] = xj;
      Itot[j][i].computedA();
    }
  }

}


void Grid::interfaceGeomFromCellPos(int j){

  double xj = Ctot[j][iLbnd[j]+1].G.x[F1];
  for (int i = 0; i < ntrack[j]-1; ++i){
    Itot[j][i].x[MV] = Ctot[j][i+1].G.x[MV] - Ctot[j][i+1].G.dx[MV]/2.;
    Itot[j][i].x[F1] = xj;
    Itot[j][i].dx[0] = Ctot[j][i].G.dx[F1];
    Itot[j][i].computedA();
  }

}


void Grid::destruct(){

  delete_array_2d(Ctot);
  delete_array_2d(Itot);

}


void Grid::apply(void (Cell::*func)()){

  // usage: grid.apply(&Cell::func_name);
  
  #pragma omp parallel for default(shared)  
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < nde_nax[MV]; ++i){
      (Ctot[j][i].*func)();
    }
  }

}


// overloading apply for FluidState
void Grid::apply(void (FluidState::*func)(), bool noborder){

  int jL, jR;
  if (noborder){ jL = 1; jR = nde_nax[F1]-1; }
  else         { jL = 0; jR = nde_nax[F1];   }

  for (int j = jL; j < jR; ++j){

    int iL, iR;
    if (noborder){ iL = 1; iR = ntrack[j]-1; }
    else         { iL = 0; iR = ntrack[j];   }

    for (int i = iL; i < iR; ++i){
      (Ctot[j][i].S.*func)();
    }
  }
    
}


void Grid::prim2cons(){

  int jL = 0; 
  int jR = nde_nax[F1];

  #pragma omp parallel for default(shared)
  for (int j = jL; j < jR; ++j){
    int iL = 0;
    int iR = ntrack[j];
    for (int i = iL; i < iR; ++i){
      double r = Ctot[j][i].G.x[r_];
      Ctot[j][i].S.prim2cons(r);
    }
  }

}


void Grid::state2flux(){

  int jL = 0; 
  int jR = nde_nax[F1];

  #pragma omp parallel for default(shared)
  for (int j = jL; j < jR; ++j){
    int iL = 0;
    int iR = ntrack[j];
    for (int i = iL; i < iR; ++i){
      double r = Ctot[j][i].G.x[r_];
      Ctot[j][i].S.state2flux(r);
    }
  }

}


void Grid::v2u(){

  // converting velocities to 4-velocities:
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]; ++i){
      double v   = 0;
      double vv[NUM_D];
      for (int d = 0; d < NUM_D; ++d){
        vv[d] = Ctot[j][i].S.prim[VV1+d];
        v += vv[d]*vv[d];
      }
      v = sqrt(v);
      double lfac = 1./sqrt(1.-v*v);
      for (int d = 0; d < NUM_D; ++d){
        Ctot[j][i].S.prim[UU1+d] *= lfac;
      }
    }
  }

}


void mpi_distribute(Grid *grid){

  int ngst = grid->ngst;
  int o = grid->origin[F1];
  int size  = grid->nde_nax[MV]-2*ngst;
  for (int j = 0; j < grid->nde_ncell[F1]; ++j){  // have to  copy track by track
    int index = o+j;
    std::copy_n(&(grid->Cinit[index][0]), size, &(grid->C[j][0]));
  }
  delete_array_2d(grid->Cinit);

  // storing index information in cell
  for (int j = 0; j < grid->nde_nax[F1]; ++j){
    for (int i = 0; i < grid->nde_nax[MV]; ++i){
      grid->Ctot[j][i].nde_id = grid->nde_nax[MV]*j + i;
      grid->Ctot[j][i].nde_ind[0] = j;
      grid->Ctot[j][i].nde_ind[1] = i;
    }
  }

}







