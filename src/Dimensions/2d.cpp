/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-16 15:14:51
*/

#include "../environment.h"
#include "../grid.h"
#include "../interface.h"
#include "../cell.h"
#include "../array_tools.h"
#include "../mpisetup.h"

static int splitGrid(int n_cell, int n_gst){
  /*
  int n_cell: total number of cells (non including ghost boundaries)
  int n_gst:  size of ghost regions.
  */
  int base = n_cell/worldsize;
  int rest = n_cell%worldsize;
  int nde_n_ax = base;
  if (worldrank<rest) nde_n_ax++;
  nde_n_ax+=2*n_gst;
  return nde_n_ax;

}

void c_grid :: initialise(s_par par){

  for (int d = 0; d < NUM_DIM; ++d) n_ax[d] = par.n_cell[d]+2*n_gst;
  nde_n_ax[MV_D_] = n_ax[MV_D_];
  nde_n_ax[FX_D1] = splitGrid(par.n_cell[FX_D1], n_gst);
  
  Ctot = array_2d<c_cell>(nde_n_ax[FX_D1],nde_n_ax[MV_D_]);
  I    = array_2d<c_interface>(nde_n_ax[FX_D1],nde_n_ax[MV_D_]+1);
  C    = &Ctot[n_gst];

}

void c_grid :: destruct(){

  delete_array_2d(C);
  delete_array_2d(I);

}

void c_grid :: print(int var){

  int size  = nde_n_ax[FX_D1] * nde_n_ax[MV_D_];
  MPI_Datatype cell_mpi = {0}; 
  generate_mpi_cell( &cell_mpi );

  if (worldrank == 0){

    c_cell **Cdump  = array_2d<c_cell>(n_ax[FX_D1], n_ax[MV_D_]);
    std::memcpy(&Cdump[0][0],&C[0][0],size * sizeof(c_cell))

    for (int j = 1; j < numprocs; ++j){
      MPI_Recv(&Cdump[index][0], size, cell_mpi, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (int j = 0; j < n_ax[FX_D1]; ++j){
      for (int i = 0; i < n_ax[MV_D_]; ++i) printf("%le \n", Cdump[j][i]);
      printf("\n");
    }
    cout << t << endl;
  }else{
    int index = noderank * nodeNy / nodesize;
    MPI_Send(&C[index][0], size, cell_mpi, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

}






