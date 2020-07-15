/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-15 16:24:45
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
    // still includes ghost cells on MV dim

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

}

void Grid::updateGhosts(){

  if (worldsize>1) { mpi_exchangeGhostTracks(); }

  // outer boundaries (OUTFLOW)
  if (worldrank==0){
    for (int j = 0; j <= jLbnd; ++j){
      ntrack[j] = ntrack[jLbnd+1];
      std::copy_n(&Ctot[jLbnd+1][0], ntrack[jLbnd+1],   &Ctot[j][0]);
      std::copy_n(&Itot[jLbnd+1][0], ntrack[jLbnd+1]-1, &Itot[j][0]);

      // udating ghost positions
      for (int i = 0; i < ntrack[j]; ++i){
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

      Itot[j][i] = Itot[j][iLbnd[j]+1];
      Itot[j][i].x[MV]   -= (iLbnd[j]-i+1) * Ctot[j][iLbnd[j]+1].G.dl[MV];
    }
    for (int i = iRbnd[j]; i < nde_nax[MV]; ++i){
      Ctot[j][i] = Ctot[j][iRbnd[j]-1];
      Ctot[j][i].G.x[MV] += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dl[MV];

      Itot[j][i-1] = Itot[j][iRbnd[j]-2];
      Itot[j][i-1].x[MV]   += (i-iRbnd[j]+1) * Ctot[j][iRbnd[j]-1].G.dl[MV];
    }
  }

}

void Grid::reconstructStates(int j, int i, int dim, int iplus, Interface *Int){

  if (Int==NULL){ Int = &Itot[j][i]; }

  #if SPATIAL_RECONSTRUCTION_ == PIECEWISE_CONSTANT_
    if (dim == MV){
      Int->SL = Ctot[j][i  ].S;
      Int->SR = Ctot[j][i+1].S;
      UNUSED(iplus);
    } else {
      Int->SL = Ctot[j][i].S;
      Int->SR = Ctot[j+1][iplus].S;
    }
  #endif

}

void Grid::computeFluxes(){

  for (int j = 0; j < nde_nax[F1]; ++j){

    // flux in MV direction (looping over interfaces)
    for (int i = 0; i < ntrack[j]-1; ++i){
      reconstructStates(j,i,MV);
      Itot[j][i].computeLambda();
      Itot[j][i].computeFlux();
    }
    for (int i = 1; i < ntrack[j]-1; ++i){
      // can't do the edges
      
      Ctot[j][i].update_dt(MV, Itot[j][i-1].lR, Itot[j][i].lL);
      for (int q = 0; q < NUM_Q; ++q){
        Ctot[j][i].flux[0][MV][q] = Itot[j][i-1].flux[q];
        Ctot[j][i].flux[1][MV][q] = Itot[j][i  ].flux[q];
      }
    }
  }
  // flux in F1 direction (building the interfaces)
  // vector of interfaces...
  for (int j = 0; j < nde_nax[F1]-1; ++j){
    int iplus = 1;
    double jpos = ( Ctot[j][1].G.x[F1] + Ctot[j+1][1].G.x[F1] )/2.;

    // resetting fluxes
    for (int i = 0; i < ntrack[j]; ++i){
      for (int q = 0; q < NUM_Q; ++q){ Ctot[j][i].flux[1][F1][q] = 0.; }
    }
    for (int i = 0; i < ntrack[j+1]; ++i){
      for (int q = 0; q < NUM_Q; ++q){ Ctot[j+1][i].flux[0][F1][q] = 0.; }
    }

    for (int i = 1; i < ntrack[j]-1; ++i){
      double xLj = Ctot[j][i].G.x[MV] - Ctot[j][i].G.dl[MV]/2.;
      double xRj = Ctot[j][i].G.x[MV] + Ctot[j][i].G.dl[MV]/2.;
      double xLjplus = -1.e15;  // unlimited ghost-cell size
      iplus--;  // looking back one cell (in case it was also contributing)
      while (xLjplus < xRj){

        double xRjplus = Itot[j+1][iplus].x[MV];

        Interface Int;
        Int.dim   = F1;
        Int.x[F1] = jpos;
        Int.x[MV] = ( fmax(xLj,xLjplus) + fmin(xRj,xRjplus) )/2.;
        Int.dl[0] = fmax(0, fmin(xRj,xRjplus) - fmax(xLj,xLjplus));
        Int.dA = Int.dl[0];

        reconstructStates(j,i,F1,iplus,&Int);
        Int.computeLambda();
        Int.computeFlux();

        Ctot[j  ][i    ].update_dt(F1, Int.lL);
        Ctot[j+1][iplus].update_dt(F1, Int.lR);
        for (int q = 0; q < NUM_Q; ++q){
          Ctot[j  ][i    ].flux[1][F1][q] += Int.flux[q];
          Ctot[j+1][iplus].flux[0][F1][q] += Int.flux[q];
        }

        xLjplus = xRjplus;
        iplus++;
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

  // do not update border cells because can lead to non-physical states
  for (int j = 1; j < nde_nax[F1]-1; ++j){
    for (int i = 1; i < ntrack[j]-1; ++i){
      Ctot[j][i].update(dt);
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

template <class T> void Grid::apply(void (T::*func)()){
    
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < nde_nax[MV]; ++i){
      (Ctot[j][i].*func)();
    }
  }

}

// overloading apply for FluidState
void Grid::apply(void (FluidState::*func)()){

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < nde_nax[MV]; ++i){
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


void mpi_distribute(Grid *grid){

  int size  = grid->nde_ncell[MV];
  for (int j = 0; j < grid->nde_ncell[F1]; ++j){  // have to  copy track by track
    int index = grid->origin[F1]+j;
    std::copy_n(&(grid->Cinit[index][0]), size, &(grid->C[j][0]));
  }
  delete_array_2d(grid->Cinit);

}







