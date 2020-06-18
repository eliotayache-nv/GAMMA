#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include <mpi.h>
#include "err.h"

#define UNUSED(x) (void)

enum{RHO,PPP,UU1,UU2,UU3};  // 3rd dimension optional so must be last
enum{DEN,TAU,SS1,SS2,SS3};
enum{x_,y_,z_};
enum{ENABLED_,DISABLED_};

// ---------------------------------------------------------------------------------------
// ENVIRONMENT OPTIONS
#define NUM_C 5               // conserved (changes with number of dimensions)
#define NUM_T 1               // tracers
#define NUM_Q (NUM_C+NUM_T)   // advected variables (tracers are placed at end of list)
#define NUM_D 2             // number of dimensions
#define MV x_            // moving dimension
#define F1 y_            // moving dimension
#define GAMMA_  (4./3.)

#define OPEN_MPI_ ENABLED_

// ---------------------------------------------------------------------------------------
// GLOBAL VARIABLES
extern int worldsize, worldrank, nodesize, noderank;
extern MPI_Comm nodecom;

// ---------------------------------------------------------------------------------------
void checkEnvironment();

struct s_par {    
  double tini;
  int n_cell[NUM_D];
};

#endif