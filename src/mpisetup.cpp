/*
* @Author: eliotayache
* @Date:   2020-06-10 15:59:03
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-22 10:39:08
*/
#include "mpi.h"
#include "err.h"
#include "environment.h"
#include "array_tools.h"
#include "fluid.h"
#include "cell.h"
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


struct s_cell
{
  int status;
  double prim[NUM_Q];
  double x[NUM_D];
  double dl[NUM_D];
};

void toStruct(Cell c, s_cell * sc)
{
  sc->status = c.status;
  arrcpy<double>(c.S.prim, sc->prim, NUM_Q);
  arrcpy<double>(c.G.x   , sc->x   , NUM_D);
  arrcpy<double>(c.G.dl  , sc->dl  , NUM_D);
}

void toClass(s_cell sc, Cell * c)
{
  c->status = sc.status;
  arrcpy<double>(sc.prim, c->S.prim, NUM_Q);
  arrcpy<double>(sc.x   , c->G.x   , NUM_D);
  arrcpy<double>(sc.dl  , c->G.dl  , NUM_D);
}

void generate_mpi_cell( MPI_Datatype * cell_mpi ){


  // // STATE MPI DATATYPE
  // c_fluid_state state;
  // MPI_Datatype state_mpi_temp, state_mpi;
  // MPI_Datatype fluxes_mpi;
  // MPI_Type_vector(NUM_D,NUM_D,NUM_Q, MPI_DOUBLE, &fluxes_mpi);
  // MPI_Type_commit(&fluxes_mpi);
  // int fcount = 3; // no. of types: only doubles
  // int fblocklengths[]={NUM_Q,NUM_Q,1};
  // MPI_Datatype ftypes[]={MPI_DOUBLE,MPI_DOUBLE,fluxes_mpi};
  // MPI_Aint foffsets[fcount];

  // foffsets[0] = (char *)&(state.prim) - (char *)(&state);
  // foffsets[1] = (char *)&(state.cons) - (char *)(&state);
  // foffsets[2] = (char *)&(state.flux) - (char *)(&state);

  // MPI_Type_create_struct(fcount,fblocklengths,foffsets,ftypes,&state_mpi_temp);
  // MPI_Type_create_resized(state_mpi_temp, foffsets[0], sizeof(state), &state_mpi);
  // MPI_Type_commit(&state_mpi);

  // // // GEOMETRY MPI DATATYPE
  // s_cell_geometry geom;
  // MPI_Datatype geom_mpi;
  // int gcount = 4; // no. of types: only doubles
  // int gblocklengths[]={1,NUM_D,NUM_D,NUM_D};
  // MPI_Datatype gtypes[]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
  // MPI_Aint goffsets[gcount];

  // goffsets[0] = (char *)&(geom.dV)  - (char *)(&geom);
  // goffsets[1] = (char *)&(geom.x)   - (char *)(&geom);
  // goffsets[2] = (char *)&(geom.dl)  - (char *)(&geom);
  // goffsets[3] = (char *)&(geom.cen) - (char *)(&geom);

  // MPI_Type_create_struct(gcount,gblocklengths,goffsets,gtypes,&geom_mpi);
  // MPI_Type_commit(&geom_mpi);

  // // CELL MPI DATATYPE
  // Cell cell;
  // int ccount = 3; // no. of types: int,state,geometry
  // int cblocklengths[]={2,1,1};
  // MPI_Datatype ctypes[]={MPI_DOUBLE,state_mpi,geom_mpi};
  // MPI_Aint coffsets[ccount];

  // coffsets[0] = (char *)&(cell.status) - (char *)(&cell);
  // coffsets[1] = (char *)&(cell.S)      - (char *)(&cell);
  // coffsets[2] = (char *)&(cell.G)      - (char *)(&cell);

  // MPI_Type_create_struct(ccount,cblocklengths,coffsets,ctypes,cell_mpi);
  // MPI_Type_commit(cell_mpi);

  // // CELL MPI DATATYPE
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

  // offsets[1] = (char *)&(sc.prim)   - (char *)(&sc);
  // offsets[2] = (char *)&(sc.x)      - (char *)(&sc);
  // offsets[3] = (char *)&(sc.dl)     - (char *)(&sc);

  printf("blah\n");

  MPI_Type_create_struct(count,blocklengths,offsets,types,cell_mpi);
  MPI_Type_commit(cell_mpi);

}



void mpi_finalise(){

  MPI_Finalize();

}
