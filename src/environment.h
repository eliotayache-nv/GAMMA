#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include <mpi.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>
#include <unistd.h>
#include <errno.h>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "err.h"

#define UNUSED(x) (void)(x)
#define ENABLED_  0
#define DISABLED_ 1

enum{RHO,PPP,UU1,UU2,UU3};  // 3rd dimension optional so must be last
enum{TP1,TP2,VV1,VV2,VV3};
enum{DEN,TAU,SS1,SS2,SS3};
enum{x_,y_,z_};      // cartesian
enum{r_,t_,p_};      // spherical
enum{X_,Y_,Z_};      // in case of typo
enum{R_,T_,P_};      // in case of typo
enum{left_,right_};    
enum{skip_,merge_,split_};

enum{PIECEWISE_CONSTANT_,PIECEWISE_LINEAR_};

// ---------------------------------------------------------------------------------------
// ENVIRONMENT OPTIONS
#define NUM_C 4               // conserved (changes with number of dimensions)
#define NUM_T 1               // tracers
#define NUM_Q (NUM_C+NUM_T)   // advected variables (tracers are placed at end of list)
#define TR1   NUM_C         // index of first tracer 
#define NUM_D 2             // number of dimensions
#define MV    x_            // moving dimension
#define F1    y_            // fixed dimension 1
#define F2    z_            // fixed dimension 2
#define VI      1.            // interface velocity (units of CD velocity)
#define GAMMA_  (5./3.)
#define CFL_    0.2

#define DUMPSTEP_ 10

#define MPI_ ENABLED_
#define OMP_ ENABLED_

#define SPATIAL_RECONSTRUCTION_ PIECEWISE_LINEAR_
#define CIRC_REGRID_ DISABLED_

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