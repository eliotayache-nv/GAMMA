#ifndef MPI_H_
#define MPI_H_

#include <mpi.h>
#include "environment.h"

extern int worldsize, worldrank, nodesize, noderank;
extern MPI_Comm nodecom;

#endif