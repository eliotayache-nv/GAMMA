#ifndef MPI_H_
#define MPI_H_

#include "environment.h"

struct s_cell;
void mpi_init(int *argc, char **argv[]);
void toStruct(Cell c, s_cell * sc);
void toClass(s_cell sc, Cell * c);
void generate_mpi_cell( MPI_Datatype * cell_mpi);
void mpi_distribute(Grid *grid);
void mpi_finalise();

#endif