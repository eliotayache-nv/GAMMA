#ifndef MPI_H_
#define MPI_H_

#include "environment.h"


class Grid;

struct s_cell
{
  int status;
  double prim[NUM_Q];
  double x[NUM_D];
  double dx[NUM_D];
  double prim0[NUM_Q];
  double x0[NUM_D];
  double dx0[NUM_D];
};

void mpi_init(int *argc, char **argv[]);
void toStruct(Cell c, s_cell * sc);
void toClass(s_cell sc, Cell * c);
void generate_mpi_cell( MPI_Datatype * cell_mpi);
void mpi_distribute(Grid *grid);
void mpi_finalise();

#endif