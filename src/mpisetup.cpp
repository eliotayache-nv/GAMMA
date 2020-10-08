/*
* @Author: eliotayache
* @Date:   2020-06-10 15:59:03
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-07 17:34:14
*/
#include "mpi.h"
#include "err.h"
#include "environment.h"
#include "array_tools.h"
#include "fluid.h"
#include "cell.h"
#include "mpisetup.h"
#include <stddef.h>
#include <iostream>
#include <cstddef>

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

void toStruct(Cell c, s_cell * sc){

  sc->status = c.status;
  arrcpy<double>(c.S.prim, sc->prim, NUM_Q);
  arrcpy<double>(c.G.x   , sc->x   , NUM_D);
  arrcpy<double>(c.G.dx  , sc->dx  , NUM_D);

}

void toClass(s_cell sc, Cell * c){

  c->status = sc.status;
  arrcpy<double>(sc.prim, c->S.prim, NUM_Q);
  arrcpy<double>(sc.x   , c->G.x   , NUM_D);
  arrcpy<double>(sc.dx  , c->G.dx  , NUM_D);
  c->computeAllGeom();
  c->S.prim2cons(c->G.x[x_]);
  c->S.state2flux(c->G.x[x_]);

}

void generate_mpi_cell( MPI_Datatype * cell_mpi ){

  s_cell sc;
  int count = 2; // no. of types: int,double
  int blocklengths[]={1,NUM_Q+NUM_D+NUM_D};
  MPI_Datatype types[]={MPI_INT,MPI_DOUBLE};
  MPI_Aint offsets[count];

  MPI_Aint a_sc, a_status, a_prim;
  MPI_Get_address(&sc,&a_sc);   // address of cell
  MPI_Get_address(&sc.status,&a_status); // address of first attribute of cell
  MPI_Get_address(&sc.prim,&a_prim); // address of first attribute of cell

  offsets[0] = a_status-a_sc;
  offsets[1] = a_prim-a_sc;

  MPI_Type_create_struct(count,blocklengths,offsets,types,cell_mpi);
  MPI_Type_commit(cell_mpi);

}




void mpi_finalise(){

  MPI_Finalize();

}
