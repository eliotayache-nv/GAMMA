/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-23 15:26:04
*/

#include "../environment.h"
#include "../grid.h"
#include "../interface.h"
#include "../cell.h"
#include "../array_tools.h"
#include "../mpisetup.h"
#include <algorithm>

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

void Grid :: initialise(s_par par){

  for (int d = 0; d < NUM_D; ++d){
    ncell[d] = par.ncell[d];
    ngst     = par.ngst;
    nax[d]   = ncell[d]+2*ngst;
  }
  nde_nax[MV]   = nax[MV];
  nde_nax[F1]   = splitGrid(ncell[F1], ngst, origin);
  nde_ncell[MV] = nde_nax[MV]-2*ngst;  // max number of active cells in mov dim
  nde_ncell[F1] = nde_nax[F1]-2*ngst;
  
  Cinit = array_2d<Cell>(ncell[F1],ncell[MV]);
  Ctot  = array_2d<Cell>(nde_nax[F1],nde_nax[MV]);
  I     = array_2d<Interface>(nde_nax[F1],nde_nax[MV]+1);
  C     = array_2d_nogst<Cell>(Ctot, nde_nax[F1], ngst);     
    // still includes ghost cells on MV dim

}

void Grid :: destruct(){

  delete_array_2d(Ctot);
  delete_array_2d(I);

}

void Grid :: print(int var){

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
      printf("coucou\n");
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
    MPI_Send(    SC,  size, cell_mpi, 0, 2, MPI_COMM_WORLD);
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







