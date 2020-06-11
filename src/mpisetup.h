#ifndef MPI_H_
#define MPI_H_

#include <mpi.h>
#include "environment.h"
#include "err.h"

extern int worldsize, worldrank, nodesize, noderank;
extern MPI_Comm nodecom;

void mpi_init(int argc, char *argv[]);

#endif