/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-26 11:24:50
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

  nsimu = 1;
  for (int d = 0; d < NUM_D; ++d){
    ncell[d] = par.ncell[d];
    ngst     = par.ngst;
    nax[d]   = ncell[d]+2*ngst;
    nsimu   *= ncell[d];
  }
  nde_nax[MV]   = nax[MV];
  nde_nax[F1]   = splitGrid(ncell[F1], ngst, origin);
  nde_ncell[MV] = nde_nax[MV]-2*ngst;  // max number of active cells in mov dim
  nde_ncell[F1] = nde_nax[F1]-2*ngst;
  nde_ntot      = nde_nax[F1] * nde_nax[MV];

  // active cells and ghost cells
  nact  = new int[nde_ncell[F1]];
  for (int j = 0; j < nde_ncell[F1]; ++j){nact[j] = nde_ncell[MV];}
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

  if (worldrank!=0 and worldrank!=worldsize-1){

    // sending to lower node, receiving from higher node
    for (int j = 0; j < ngst; ++j){
      int jout = j+ngst;
      int jin  = jRbnd+j;
      int nout = nact[jout]+2*ngst; // number of cells to receive
      int nin;                        // number of cells to send
      int tag1 = j;
      int tag2 = ngst+j;
      MPI_Sendrecv(&nout, 1, MPI_INT, rL, tag1, 
                   &nin , 1, MPI_INT, rR, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      s_cell *SCout = new s_cell[nout];
      s_cell *SCin  = new s_cell[nin];
      for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
      MPI_Sendrecv(SCout, nout, cell_mpi, rL, tag2, 
                   SCin , nin , cell_mpi, rR, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }    
    }

    // sending to higher node, receiving from lower node
    for (int j = 0; j < ngst; ++j){
      int jout = jRbnd-ngst+j;
      int jin  = j;
      int nout = nact[jout]+2*ngst; // number of cells to receive
      int nin;                        // number of cells to send
      int tag1 = j;
      int tag2 = ngst+j;
      MPI_Sendrecv(&nout, 1, MPI_INT, rR, tag1, 
                   &nin , 1, MPI_INT, rL, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      s_cell *SCout = new s_cell[nout];
      s_cell *SCin  = new s_cell[nin];
      for (int i = 0; i < nout; ++i) { toStruct(Ctot[jout][i], &SCout[i]); }
      MPI_Sendrecv(SCout, nout, cell_mpi, rR, tag2, 
                   SCin , nin , cell_mpi, rL, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i < nin; ++i) { toClass(SCin[i], &Ctot[jin][i]); }    
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

}

void Grid::updateGhosts(){

  mpi_exchangeGhostTracks();

  // outer boundaries (OUTFLOW)
  if (worldrank==0){
    for (int j = 0; j <= jLbnd; ++j){
      std::copy_n(&C[0][0], nact[j], &Ctot[j][0]);
    }
  }
  if (worldrank==worldsize-1){
    for (int j = jRbnd; j < nde_nax[F1]; ++j){
      std::copy_n(&C[nde_ncell[F1]-1][0], nact[j], &Ctot[j][0]);
    }
  }
  for (int j = 0; j < nde_nax[F1]; ++j)
  {
    for (int i = 0; i <= iLbnd[j]; ++i){Ctot[j][i] = Ctot[j][iLbnd[j]+1];}
    for (int i = iRbnd[j]; i < nde_nax[MV]; ++i){Ctot[j][i] = Ctot[j][iRbnd[j]-1];}
  }

}

void Grid::interfaceGeomFromCellPos(){

  for (int j = 0; j < nde_nax[F1]; ++j){
    double xj = C[0][0].G.x[F1];
    for (int i = 0; i < nde_nax[MV]-1; ++i){
      Itot[j][i].x[MV] = (Ctot[j][i].G.x[MV] + Ctot[j][i+1].G.x[MV])/2.;
      Itot[j][i].x[F1] = xj;
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

template <> void Grid::apply<FluidState>(void (FluidState::*func)()){

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
        printf("%le \n", Cdump[j][i].S.prim[var]);
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







