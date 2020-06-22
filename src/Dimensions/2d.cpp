/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-22 16:02:43
*/

#include "../environment.h"
#include "../grid.h"
#include "../interface.h"
#include "../cell.h"
#include "../array_tools.h"
#include "../mpisetup.h"
#include <algorithm>

static int splitGrid(int n_cell, int n_gst, int *origin){
  /*
  int n_cell:   total number of cells (non including ghost boundaries)
  int n_gst:    size of ghost regions.
  int *origin:  coordinates of C[0][0] in the simulation domain
  */
  int base = n_cell/worldsize;
  int rest = n_cell%worldsize;
  int nde_n_ax = base;
  if (worldrank<rest) nde_n_ax++;
  nde_n_ax+=2*n_gst;

  origin[MV] = 0;
  origin[F1] = worldrank*base;
  if (worldrank <rest) origin[F1]+=worldrank;
  else origin[F1]+=rest;

  return nde_n_ax;

}

void Grid :: initialise(s_par par){

  for (int d = 0; d < NUM_D; ++d){
    n_cell[d] = par.n_cell[d];
    n_gst     = par.n_gst;
    n_ax[d]   = n_cell[d]+2*n_gst;
  }
  nde_n_ax[MV]   = n_ax[MV];
  nde_n_ax[F1]   = splitGrid(n_cell[F1], n_gst, origin);
  nde_n_cell[MV] = nde_n_ax[MV]-2*n_gst;  // max number of active cells in mov dim
  nde_n_cell[F1] = nde_n_ax[F1]-2*n_gst;
  
  Cinit = array_2d<Cell>(n_cell[F1],n_cell[MV]);
  Ctot  = array_2d<Cell>(nde_n_ax[F1],nde_n_ax[MV]);
  I     = array_2d<Interface>(nde_n_ax[F1],nde_n_ax[MV]+1);
  C     = &Ctot[n_gst];     // still includes ghost cells on MV dim

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
    Cell    **Cdump = array_2d<Cell>(n_cell[F1], n_ax[MV]);
    s_cell **SCdump = array_2d<s_cell>(n_cell[F1], n_ax[MV]);

    sizes[0] = nde_n_cell[F1] * nde_n_ax[MV];
    std::copy_n(&C[0][0], sizes[0], &Cdump[0][0]);

    for (int j = 1; j < worldsize; ++j){
      int o[NUM_D];
      MPI_Recv(            &sizes[j],        1,  MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(                    o,    NUM_D,  MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&SCdump[o[F1]][o[MV]], sizes[j], cell_mpi, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (int j = 0; j < n_cell[F1]; ++j){
      for (int i = 0; i < n_ax[MV]; ++i) {
        toClass(SCdump[j][i], &Cdump[j][i]);
        printf("%le \n", Cdump[j][i].S.prim[var]);
      }
      printf("\n");
    }
  }else{
    int size  = nde_n_cell[F1] * nde_n_ax[MV];  // size includes MV ghost cells
    MPI_Send( &size,     1,  MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(origin, NUM_D,  MPI_INT, 0, 1, MPI_COMM_WORLD);
    s_cell **SC = array_2d<s_cell>(nde_n_cell[F1],nde_n_ax[MV]);
    for (int j = 0; j < nde_n_cell[F1]; ++j){
      for (int i = 0; i < nde_n_ax[MV]; ++i){
        toStruct(C[j][i], &SC[j][i]);
      }
    }
    MPI_Send(    SC,  size, cell_mpi, 0, 2, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

}


void mpi_distribute(Grid *grid){

  int size  = grid->nde_n_ax[MV];
  for (int j = 0; j < grid->nde_n_cell[F1]; ++j){  // have to  copy track by track
    int index = grid->origin[F1]+j;
    std::copy_n(&(grid->Cinit[index]), size, &(grid->C[j]));
  }
  delete_array_2d(grid->Cinit);

}







