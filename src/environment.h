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
#define PIECEWISE_CONSTANT_ 0
#define PIECEWISE_LINEAR_   1
#define IDEAL_EOS_      0
#define SYNGE_EOS_      1

enum{RHO,PPP,UU1,UU2,UU3};  // 3rd dimension optional so must be last
enum{TP1,TP2,VV1,VV2,VV3};
enum{DEN,TAU,SS1,SS2,SS3};
enum{x_,y_,z_};      // cartesian
enum{r_,t_,p_};      // spherical
enum{X_,Y_,Z_};      // in case of typo
enum{R_,T_,P_};      // in case of typo
enum{left_,right_};    
enum{skip_,merge_,split_};


// ---------------------------------------------------------------------------------------
// ENVIRONMENT OPTIONS
#define NUM_C  4             // conserved (changes with number of dimensions)
#define NUM_TR 2             // user-specified tracers
#define NUM_D  1             // number of dimensions
#define MV     x_            // moving dimension
#define F1     y_            // fixed dimension 1
#define F2     z_            // fixed dimension 2
#define VI     1.            // interface velocity (units of CD velocity)
#define GAMMA_ (5./3.)
#define EOS_   SYNGE_EOS_
#define CFL_   0.2

#define DUMPSTEP_ 100

#define MPI_ ENABLED_
#define OMP_ ENABLED_

#define SPATIAL_RECONSTRUCTION_ PIECEWISE_LINEAR_
#define CIRC_REGRID_            ENABLED_
#define SHOCK_DETECTION_        ENABLED_
#define DETECT_SHOCK_THRESHOLD_ 0.1
#define LOCAL_SYNCHROTRON_      ENABLED_
#define GAMMA_MAX_INIT_         (1.e7)

// ---------------------------------------------------------------------------------------
// DO NOT MODIFY!
#if LOCAL_SYNCHROTRON_ == ENABLED_
  #define NUM_S 2
  #define GMN   (NUM_C+NUM_TR)    // gamma_min
  #define GMX   (NUM_C+NUM_TR+1)  // gamma_max
#else
  #define NUM_S 0
#endif 

#define TR1   NUM_C                 // index of first tracer 
#define NUM_T (NUM_TR+NUM_S)        // total num of tracers
#define NUM_Q (NUM_C+NUM_T)   // advected variables (tracers are placed at end of list)

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