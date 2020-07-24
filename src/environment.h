#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include <mpi.h>
#include <math.h>
#include <vector>
#include "err.h"

#define UNUSED(x) (void)(x)

enum{RHO,PPP,UU1,UU2,UU3};  // 3rd dimension optional so must be last
enum{TP1,TP2,VV1,VV2,VV3};
enum{DEN,TAU,SS1,SS2,SS3};
enum{x_,y_,z_};
enum{X_,Y_,Z_}; // in case of typo
enum{ENABLED_,DISABLED_};

enum{PIECEWISE_CONSTANT_,PIECEWISE_LINEAR_};

// ---------------------------------------------------------------------------------------
// ENVIRONMENT OPTIONS
#define NUM_C 4               // conserved (changes with number of dimensions)
#define NUM_T 1               // tracers
#define NUM_Q (NUM_C+NUM_T)   // advected variables (tracers are placed at end of list)
#define NUM_D 2             // number of dimensions
#define MV    x_            // moving dimension
#define F1    y_            // flixed dimension
#define VI    1.            // interface velocity (units of CD velocity)
#define GAMMA_  (5./3.)
#define CFL_    0.2

#define DUMPSTEP_ 10

#define OPEN_MPI_ ENABLED_

#define SPATIAL_RECONSTRUCTION_ PIECEWISE_CONSTANT_

// ---------------------------------------------------------------------------------------
// GLOBAL VARIABLES
extern int worldsize, worldrank, nodesize, noderank;
extern MPI_Comm nodecom;

// ---------------------------------------------------------------------------------------
void checkEnvironment();

struct s_par {    
  double tini;
  int ncell[NUM_D];
  int ngst;
  int nmax;
};

#endif