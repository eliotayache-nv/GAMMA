#ifndef MPI_H_
#define MPI_H_

#include "environment.h"

void mpi_init(int *argc, char **argv[]);
void generate_mpi_cell( MPI_Datatype * cell_mpi);
void mpi_finalise();

#endif