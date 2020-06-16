/*
* @Author: eliotayache
* @Date:   2020-06-10 15:59:03
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-16 14:58:49
*/
#include "mpi.h"
#include "err.h"
#include "environment.h"
#include "fluid.h"
#include "cell.h"
#include <stddef.h>

void mpi_init(int *argc, char **argv[]){

  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
  MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up

  // get node-specific information
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, worldrank,
    MPI_INFO_NULL, &nodecom);

  MPI_Comm_size(nodecom, &nodesize);
  MPI_Comm_rank(nodecom, &noderank);

  // checking that the number of tasks is equal to the number of nodes
  // We should have exactly one process per node
  // That way we can use openMP on each node and not have to worry about splitting
  // the grid manually for each process
  int nodesize_check=1;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&nodesize, &nodesize_check, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    // to be able to throw exception on all processes simultaneously
  if (nodesize_check!=1) throw MPITooManyTasksException();

}

void generate_mpi_cell( MPI_Datatype * cell_mpi ){

  // STATE MPI DATATYPE
  c_fluid_state state;
  MPI_Datatype state_mpi;
  int fcount = 1; // no. of types: only doubles
  int fblocklenghts[]={NUM_Q*NUM_Q*NUM_DIM*NUM_Q};
  MPI_Datatype ftypes[]={MPI_DOUBLE};
  MPI_Aint foffsets[fcount];

  foffsets[0] = (char *)&(state.prim) - (char *)(&state);

  MPI_Type_create_struct(fcount,fblocklenghts,foffsets,ftypes,&state_mpi);
  MPI_Type_commit(&state_mpi);

  // GEOMETRY MPI DATATYPE
  MPI_Datatype geom_mpi;
  int gcount = 1; // no. of types: only doubles
  int gblocklenghts[]={1+4*NUM_DIM};
  MPI_Datatype gtypes[]={MPI_DOUBLE};
  MPI_Aint goffsets[gcount];

  goffsets[0] = offsetof(s_cell_geometry,dV);

  MPI_Type_create_struct(gcount,gblocklenghts,goffsets,gtypes,&geom_mpi);
  MPI_Type_commit(&geom_mpi);

  // CELL MPI DATATYPE
  c_cell cell;
  int ccount = 3; // no. of types: int,state,geometry
  int cblocklenghts[]={2,1,1};
  MPI_Datatype ctypes[]={MPI_DOUBLE,state_mpi,geom_mpi};
  MPI_Aint coffsets[ccount];

  coffsets[0] = (char *)&(cell.status) - (char *)(&cell);
  coffsets[0] = (char *)&(cell.S)      - (char *)(&cell);
  coffsets[0] = (char *)&(cell.G)      - (char *)(&cell);

  MPI_Type_create_struct(ccount,cblocklenghts,coffsets,ctypes,cell_mpi);
  MPI_Type_commit(cell_mpi);
}


void mpi_finalise(){

  MPI_Finalize();

}
