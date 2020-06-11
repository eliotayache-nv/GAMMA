/*
* @Author: eliotayache
* @Date:   2020-06-10 15:59:03
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 18:34:57
*/
#include "mpi.h"
#include "err.h"
#include "environment.h"

void mpi_init(int *argc, char **argv[])
{
  #if OPEN_MPI_ == ENABLED_

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

  #elif OPEN_MPI_ == DISABLED_
    worldsize = 1;
    worldrank = 0;
    noderank  = 0;
  #endif   
}

void mpi_finalise()
{
  MPI_Finalize();
}
