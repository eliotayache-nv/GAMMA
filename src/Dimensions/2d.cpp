/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-08-26 09:17:47
*/

#include "../environment.h"
#include "../grid.h"
#include "../interface.h"
#include "../cell.h"
#include "../array_tools.h"
#include "../mpisetup.h"
#include <algorithm>


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
  
  Cinit = array_2d<Cell>(ncell[F1],ncell[MV]);
  Ctot  = array_2d<Cell>(nde_nax[F1],nde_nax[MV]);
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
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int j = 0; j < ngst; ++j){
    for (int i = 0; i < nde_nax[MV]; ++i){

      // updating interface positions from communicated cells
      if (i != nde_nax[MV]-1){ 
        int jn = nde_nax[F1]-1-j;
        Itot[j][i].x[MV]  = Ctot[j][i].G.x[MV] + Ctot[j][i].G.dl[MV]/2.;
        Itot[jn][i].x[MV] = Ctot[jn][i].G.x[MV] + Ctot[jn][i].G.dl[MV]/2.;
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


void Grid::updateGhosts(){

  if (worldsize>1) { mpi_exchangeGhostTracks(); }

  // outer boundaries (OUTFLOW)
  if (worldrank==0){
    for (int j = 0; j <= jLbnd; ++j){
      ntrack[j] = ntrack[jLbnd+1];
      std::copy_n(&Ctot[jLbnd+1][0], ntrack[jLbnd+1],   &Ctot[j][0]);
      std::copy_n(&Itot[jLbnd+1][0], ntrack[jLbnd+1]-1, &Itot[j][0]);

      // udating ghost positions, ids and indexes
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] -= (jLbnd-j+1)*C[0][0].G.dl[F1];
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] -= (jLbnd-j+1)*C[0][0].G.dl[F1]; 
        }
      }
    }
  }
  if (worldrank==worldsize-1){
    for (int j = jRbnd; j < nde_nax[F1]; ++j){
      ntrack[j] = ntrack[jRbnd-1];
      std::copy_n(&Ctot[jRbnd-1][0], ntrack[jRbnd-1]  , &Ctot[j][0]);
      std::copy_n(&Itot[jRbnd-1][0], ntrack[jRbnd-1]-1, &Itot[j][0]);

      // udating ghost positions
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dl[F1];
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dl[F1];
        }
      }
    }
  }
  for (int j = 0; j < nde_nax[F1]; ++j)
  {
    for (int i = 0; i <= iLbnd[j]; ++i){
      Ctot[j][i] = Ctot[j][iLbnd[j]+1];
      Ctot[j][i].G.x[MV] -= (iLbnd[j]-i+1) * Ctot[j][iLbnd[j]+1].G.dl[MV];
      int ind[] = {j,i};
      assignId(ind);

      Itot[j][i] = Itot[j][iLbnd[j]+1];
      Itot[j][i].x[MV]   -= (iLbnd[j]-i+1) * Ctot[j][iLbnd[j]+1].G.dl[MV];
    }
    for (int i = iRbnd[j]; i < nde_nax[MV]; ++i){
      Ctot[j][i] = Ctot[j][iRbnd[j]-1];
      Ctot[j][i].G.x[MV] += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dl[MV];
      int ind[] = {j,i};
      assignId(ind);

      Itot[j][i-1] = Itot[j][iRbnd[j]-2];
      Itot[j][i-1].x[MV] += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dl[MV];
    }
  }

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
      grad[q] = (qR - qL) / (xR -xL);
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
      dln_mv = c_out->G.dl[MV];
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
          IS[s]->prim2cons();
          IS[s]->state2flux();
        }
      }
      UNUSED(idn);

    } else {

      if (j == 0  or j == nde_nax[F1]-2){
        Int->SL = Ctot[j][i].S;
        Int->SR = Ctot[0][idn].S;
      }
      else {
        // recover gardients in moving direction and compute states at projected location
        // careful with the boundaries for LL and RR
        Cell *cL = &Ctot[j][i]; 
        Cell *cR = &Ctot[0][idn];
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
          Int->SL.prim2cons();
          Int->SR.prim2cons();
          Int->SL.state2flux();
          Int->SR.state2flux();
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
        c->neigh[d][1].clear();

        if (d == MV){
          c->neigh[d][0].push_back(Ctot[j][i-1].nde_id);
          c->neigh[d][1].push_back(Ctot[j][i+1].nde_id);
        } else {
          double xjL = c->G.x[MV] - c->G.dl[MV]/2.;
          double xjR = c->G.x[MV] + c->G.dl[MV]/2.;

          if (j!=0){
            double xm = Itot[j-1][im].x[MV];
            if (xm > xjL and im > 0){ 
              im--; xm = Itot[j-1][im].x[MV];
            }
            while (xm < xjR){
              c->neigh[d][0].push_back(Ctot[j-1][im+1].nde_id);
              if (im > ntrack[j]-1) break;
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
              if (ip > ntrack[j]-1) break;
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
  

void Grid::movDir_ComputeLambda(){

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      reconstructStates(j,i,MV);
      Itot[j][i].computeLambda();
    }
  }

}

void Grid::updateKinematics(){

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      double v = VI * Itot[j][i].lS;
      double lfac = 1./sqrt(1.- v*v);
      Itot[j][i].v = v;
      Itot[j][i].lfac = lfac;
    }
  }  

}

void Grid::computeFluxes(){

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
      double xL0 = x0 - c0->G.dl[MV]/2.;
      double xR0 = x0 + c0->G.dl[MV]/2.;

      for (std::vector<int>::size_type n = 0; n < c0->neigh[F1][1].size(); ++n){
        int idn = c0->neigh[F1][1][n];
        Cell *cn = &Ctot[0][idn];
        double xn = cn->G.x[MV];
        double xLn = xn - cn->G.dl[MV]/2.;
        double xRn = xn + cn->G.dl[MV]/2.;

        Interface Int;
        Int.dim   = F1;
        Int.x[F1] = jpos;
        Int.x[MV] = ( fmax(xL0,xLn) + fmin(xR0,xRn) )/2.;
        Int.dl[0] = fmax(0, fmin(xR0,xRn) - fmax(xL0,xLn));
        Int.dA = Int.dl[0];

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
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]; ++i){
      dt = fmin(dt, Ctot[j][i].dt_loc);
      Ctot[j][i].dt_loc = 1.e15;  // resetting dt_loc for next update
    }
  }

  return dt;

}


void Grid::update(double dt){

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ntrack[j]-1; ++i){
      Itot[j][i].move(dt);
    }
  }
  // do not update border cells because can lead to non-physical states
  for (int j = 1; j < nde_nax[F1]-1; ++j){
    for (int i = 1; i < ntrack[j]-1; ++i){
      double xL = Itot[j][i-1].x[MV];
      double xR = Itot[j][i].x[MV];
      // printf("%d %d\n", j, i);
      Ctot[j][i].update(dt,xL,xR);
    }
  }

}


void Grid::interfaceGeomFromCellPos(){

  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    double xj = Ctot[j][iLbnd[j]+1].G.x[F1];
    for (int i = 0; i < ntrack[j]-1; ++i){
      Itot[j][i].x[MV] = (Ctot[j][i].G.x[MV] + Ctot[j][i+1].G.x[MV])/2.;
      Itot[j][i].x[F1] = xj;
      Itot[j][i].dl[0] = Ctot[j][i].G.dl[F1];
      Itot[j][i].dA = Itot[j][i].dl[0];
    }
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
      // printf("%d %d\n", j, i);
      (Ctot[j][i].S.*func)();
    }
  }
    
}

void Grid::print(int var){

  MPI_Datatype cell_mpi = {0}; 
  generate_mpi_cell(&cell_mpi);

  if (worldrank == 0){

    int sizes[worldsize];
    Cell    **Cdump = array_2d<Cell>(nax[F1], nax[MV]);
    s_cell **SCdump = array_2d<s_cell>(nax[F1], nax[MV]);
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i < nde_nax[MV]; ++i){
        toStruct(Ctot[j][i], &SCdump[j][i]);
      }
    }

    sizes[0] = nde_nax[F1] * nde_nax[MV];
    std::copy_n(&Ctot[0][0], sizes[0], &Cdump[0][0]);

    for (int j = 1; j < worldsize; ++j){
      int o[NUM_D]; // origin
      MPI_Recv(      &sizes[j],        1,  MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(              o,    NUM_D,  MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int i1 = o[F1];
      int i2 = o[MV];
      MPI_Recv(&SCdump[i1][i2], sizes[j], cell_mpi, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (int j = ngst; j < ncell[F1]+ngst; ++j){
      for (int i = ngst; i < ncell[MV]+ngst; ++i) {
        toClass(SCdump[j][i], &Cdump[j][i]);
        printf("%le ", Cdump[j][i].S.prim[var]);
      }
      printf("\n");
    }

  }else{
    int size  = nde_nax[F1] * nde_nax[MV];  // size includes MV ghost cells
    MPI_Send( &size,     1,  MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(origin, NUM_D,  MPI_INT, 0, 1, MPI_COMM_WORLD);
    s_cell **SC = array_2d<s_cell>(nde_nax[F1],nde_nax[MV]);
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i < nde_nax[MV]; ++i){
        toStruct(Ctot[j][i], &SC[j][i]);
      }
    }
    MPI_Send(&SC[0][0],  size, cell_mpi, 0, 2, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

}


void Grid::printCols(int var){

  MPI_Datatype cell_mpi = {0}; 
  generate_mpi_cell(&cell_mpi);

  if (worldrank == 0){

    int sizes[worldsize];
    Cell    **Cdump = array_2d<Cell>(nax[F1], nax[MV]);
    s_cell **SCdump = array_2d<s_cell>(nax[F1], nax[MV]);
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i < nde_nax[MV]; ++i){
        toStruct(Ctot[j][i], &SCdump[j][i]);
      }
    }

    sizes[0] = nde_nax[F1] * nde_nax[MV];
    std::copy_n(&Ctot[0][0], sizes[0], &Cdump[0][0]);

    for (int j = 1; j < worldsize; ++j){
      int o[NUM_D]; // origin
      MPI_Recv(      &sizes[j],        1,  MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(              o,    NUM_D,  MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int i1 = o[F1];
      int i2 = o[MV];
      MPI_Recv(&SCdump[i1][i2], sizes[j], cell_mpi, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("x, y, z\n");
    for (int j = ngst; j < ncell[F1]+ngst; ++j){
      for (int i = ngst; i < ncell[MV]+ngst; ++i) {
        toClass(SCdump[j][i], &Cdump[j][i]);
        printf("%le %le %le\n", 
          Cdump[j][i].G.x[x_],
          Cdump[j][i].G.x[y_],
          Cdump[j][i].S.prim[var]);
      }
    }

  }else{
    int size  = nde_nax[F1] * nde_nax[MV];  // size includes MV ghost cells
    MPI_Send( &size,     1,  MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(origin, NUM_D,  MPI_INT, 0, 1, MPI_COMM_WORLD);
    s_cell **SC = array_2d<s_cell>(nde_nax[F1],nde_nax[MV]);
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i < nde_nax[MV]; ++i){
        toStruct(Ctot[j][i], &SC[j][i]);
      }
    }
    MPI_Send(&SC[0][0],  size, cell_mpi, 0, 2, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

}


void mpi_distribute(Grid *grid){

  int size  = grid->nde_ncell[MV];
  for (int j = 0; j < grid->nde_ncell[F1]; ++j){  // have to  copy track by track
    int index = grid->origin[F1]+j;
    std::copy_n(&(grid->Cinit[index][0]), size, &(grid->C[j][0]));
  }
  delete_array_2d(grid->Cinit);

  // storing index information in cell
  for (int j = 0; j < grid->nde_nax[F1]; ++j){
    for (int i = 0; i < grid->nde_nax[MV]; ++i){
      grid->Ctot[j][i].nde_id = grid->nde_nax[MV]*j + i;
      grid->Ctot[j][i].nde_ind[0] = j;
      grid->Ctot[j][i].nde_ind[1] = i;
      // if(j==101) printf("%d %d %d\n",
      //   grid->Ctot[j][i].nde_id,
      //   grid->Ctot[j][i].nde_ind[0] = j,
      //   grid->Ctot[j][i].nde_ind[1] = i);
    }
  }

}







