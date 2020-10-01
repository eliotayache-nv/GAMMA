/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-01 18:25:51
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
    delete(SCin); delete(SCout);
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
    delete(SCin); delete(SCout);
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
        }
      }
    }
  }
  if (worldrank==worldsize-1){
    for (int j = jRbnd; j < nde_nax[F1]; ++j){
      ntrack[j] = ntrack[jRbnd-1];
      nact[j] = nact[jRbnd-1];
      iRbnd[j] = iRbnd[jRbnd-1];
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

      Itot[j][i] = Itot[j][iLbnd[j]+1];
      Itot[j][i].x[MV]   -= (iLbnd[j]-i+1) * Ctot[j][iLbnd[j]+1].G.dx[MV];
    }
    for (int i = iRbnd[j]; i < nde_nax[MV]; ++i){
      Ctot[j][i] = Ctot[j][iRbnd[j]-1];
      Ctot[j][i].G.x[MV] += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dx[MV];
      Ctot[j][i].computeAllGeom();
      int ind[] = {j,i};
      assignId(ind);

      Itot[j][i-1] = Itot[j][iRbnd[j]-2];
      Itot[j][i-1].x[MV] += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dx[MV];
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

    for (int q = 0; q < NUM_Q; ++q){
      // reconstruction using primitive variables
      double qL = cL.S.prim[q];
      double qR = cR.S.prim[q];
      double xL = cL.G.cen[dim];
      double xR = cR.G.cen[dim];
      grad[q] = (qR - qL) / (xR-xL);
    }

  }

  static Cell* findNeighbor(Cell *c, int side, double x_mv, Grid *g){

    Cell *c_out;
    Cell **Ctot = g->Ctot;
    int n = 0;
    int idn;
    int n_neigh = c->neigh[F1][side].size();
    double xn_mv, dln_mv;
    do {
      idn = c->neigh[F1][side][n];
      c_out  = &Ctot[0][idn];
      xn_mv  = c_out->G.x[MV];
      dln_mv = c_out->G.dx[MV];
      n++;
    } while (xn_mv + 0.5*dln_mv < x_mv and n < n_neigh);
      // if no aligned neighbor, we return state from closest left neighbor to the left

    return c_out;

  }

#endif


// void Grid::gradients(Cell *c){

//   for (int d = 0; d < NUM_D; ++d){
//     int indL[NUM_D];
//     int indR[NUM_D];
//     indL[d] = c->nde_ind[d] - 1;
//     indR[d] = c->nde_ind[d] + 1;
//     Cell cL = Ctot[indL[F1]][indL[MV]];
//     Cell cR = Ctot[indR[F1]][indR[MV]];
//     grad(cL, cR, d, &(c->grad));
//   }

// }


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
            double xc = c[s]->G.x[dim];
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
          or cR->neigh[F1][0].size() == 0){
        Int->SL = Ctot[j][i].S;
        Int->SR = Ctot[0][idn].S;
      }
      else {
        // recover gardients in moving direction and compute states at projected location
        // careful with the boundaries for LL and RR
        Cell *cLL, *cRR;
        double xm = Int->x[MV]; 
        double xf = Int->x[F1];
        cLL = findNeighbor(cL, left_,  xm, this);
        cRR = findNeighbor(cR, right_, xm, this);
        double dxmL  = xm - cL->G.x[MV];
        double dxmR  = xm - cR->G.x[MV];
        double dxmLL = xm - cLL->G.x[MV];
        double dxmRR = xm - cRR->G.x[MV];
        double dxfL  = xf - cL->G.x[F1];
        double dxfR  = xf - cR->G.x[F1];

        for (int q = 0; q < NUM_Q; ++q){
          double qL = cL->S.prim[q] + cL->grad_mv[q] * dxmL;
          double qR = cR->S.prim[q] + cR->grad_mv[q] * dxmR;
          double qLL = cLL->S.prim[q] + cLL->grad_mv[q] * dxmLL;
          double qRR = cRR->S.prim[q] + cRR->grad_mv[q] * dxmRR;
          double gradL = (qL - qLL) / (cL->G.x[F1] - cLL->G.x[F1]);
          double grad0 = (qR - qL)  / (cR->G.x[F1] - cL->G.x[F1]);
          double gradR = (qR - qRR) / (cR->G.x[F1] - cRR->G.x[F1]);
          gradL = minmod(gradL, grad0);
          gradR = minmod(grad0, gradR);
          Int->SL.prim[q] = qL + gradL * dxfL;
          Int->SR.prim[q] = qR + gradR * dxfR;
          Int->SL.prim2cons(Int->x[x_]);
          Int->SR.prim2cons(Int->x[x_]);
          Int->SL.state2flux(Int->x[x_]);
          Int->SR.state2flux(Int->x[x_]);
        }
      }
    }
    // double gradL[NUM_Q], gradR[NUM_Q];
    // if (dim == MV){
    //   Cell *cL  = &Ctot[j][i  ];
    //   Cell *cR  = &Ctot[j][i+1];
    //   gradients(cL); gradients(cR);
    //   Int->SL = cL->S;
    //   Int->SR = cR->S;
    //   UNUSED(iplus);
    // } else {
    //   Int->SL = Ctot[j][i].S;
    //   Int->SR = Ctot[j+1][iplus].S;
    // }
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
              if (im > ntrack[j-1]-1) break;
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
              ip++; 
              if (ip > ntrack[j+1]-1) break;
              xp = Itot[j+1][ip].x[MV];
            }
          }
        }
      }
    }
  }

  if (print){
    for (int j = 1; j < nde_nax[F1]-1; ++j){
      for (int i = 1; i < ntrack[j]-1; ++i){
        for (int d = 0; d < NUM_D; ++d){
          for (int s = 0; s < 2; ++s)
          {
            printf("%d %d, %d %d: ", j, i, d, s);
            Cell c  = Ctot[j][i];

            for (std::vector<int>::size_type n = 0; n < c.neigh[d][s].size(); ++n){

              int id  = c.neigh[d][s][n];
              Cell cn = Ctot[0][id];
              int jn  = cn.nde_ind[0];
              int in  = cn.nde_ind[1];
              printf(" %d %d", jn, in);
            }
            printf("\n");  
          }
        }
      }
    } 
  }

}


void Grid::regrid(){

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    targetRegridVictims(j); // updating ismall and ibig
    for (int i = iLbnd[j]+2; i <= iRbnd[j]-2; ++i){  // only active cells
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
  for (int i = iLbnd[j]+2; i <= iRbnd[j]-2; ++i){   // not allowing edges to be victims
    Cell c = Ctot[j][i];

    double val = c.regridVal();
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

  // Simple copy so far
  for (int q = 0; q < NUM_Q; ++q){
    cR->S.prim[q] = c->S.prim[q];
  }

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
  if (cL->G.dx[MV] < cR->G.dx[MV]){
    side = left_;
    cVic = cL;
  } else {
    side = right_;
    cVic = cR;
  }

  // updating values
  double dV_loc = c->G.dV;
  double dV_vic = cVic->G.dV;
  for (int q = 0; q < NUM_Q; ++q){
    double Q_loc = c->S.prim[q];
    double Q_vic = cVic->S.prim[q];
    c->S.prim[q] = (Q_loc * dV_loc + Q_vic * dV_vic) / (dV_loc + dV_vic);
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

  c->S.prim2cons(c->G.x[r_]);
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
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      reconstructStates(j,i,MV);
      Itot[j][i].computeLambda();
    }
  }

}

void Grid::updateKinematics(){

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      double v = VI * Itot[j][i].lS;
      double lfac = 1./sqrt(1.- v*v);
      Itot[j][i].v = v;
      Itot[j][i].lfac = lfac;
    }
  }  
  userKinematics();   // overriding with user preferences

}

void Grid::computeFluxes(){

  #pragma omp parallel for default(shared)
  for (int j = 0; j < nde_nax[F1]; ++j){

    // flux in MV direction (looping over interfaces)
    for (int i = 0; i < ntrack[j]-1; ++i){
      Itot[j][i].computeFlux();
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
    double jpos = ( Ctot[j][1].G.x[F1] + Ctot[j+1][1].G.x[F1] )/2.;


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
  for (int j = 0; j < nde_nax[F1]; ++j){
    // int proc = omp_get_thread_num();
    // int nthrd = omp_get_num_threads();
    for (int i = 0; i < ntrack[j]-1; ++i){
      // printf("%d from %d of %d on node %d of %d\n", j,proc, nthrd, worldrank, worldsize);
      Itot[j][i].move(dt);
      // if (i==ntrack[j]-2) printf("%le\n", Itot[j][i].x[MV]);
    }
  }
  // exit(10);
  // do not update border cells because can lead to non-physical states
  #pragma omp parallel for default(shared)
  for (int j = 1; j < nde_nax[F1]-1; ++j){
    // int proc = omp_get_thread_num();
    // int nthread = omp_get_num_threads();
    // printf("%d hello from %d of %d on node %d from %d\n",j, proc, nthread, worldrank, worldsize );
    for (int i = 1; i < ntrack[j]-1; ++i){
      double xL = Itot[j][i-1].x[MV];
      double xR = Itot[j][i].x[MV];
      Ctot[j][i].update(dt,xL,xR);
    }
  }

}


void Grid::interfaceGeomFromCellPos(){

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    double xj = Ctot[j][iLbnd[j]+1].G.x[F1];
    for (int i = 0; i < ntrack[j]-1; ++i){
      Itot[j][i].x[MV] = Ctot[j][i+1].G.x[MV] - Ctot[j][i+1].G.dx[MV]/2.;
      Itot[j][i].x[F1] = xj;
      Itot[j][i].dx[0] = Ctot[j][i].G.dx[F1];
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







