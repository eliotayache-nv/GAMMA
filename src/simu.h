#ifndef SIM_H_
#define SIM_H_

#include "environment.h"
#include "grid.h"

void loadParams(s_par *par);

class Simu
{
public:
    Simu();
    ~Simu();
    
    bool      stop;   // stop computation
    double    t, dt, tmax;    // time, timestep, max simulation time
    long int  it, itmax;      // time increment, max time increment
    Grid      grid;           // simulation grid
    s_par     par;

    FILE*     fout;

    int initialise();   // load configurations
    int run();          // Runs the hydro evolution of the grid 

};

#endif