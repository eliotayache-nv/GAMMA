#ifndef SIM_H_
#define SIM_H_

#include "environment.h"
#include "grid.h"

struct s_par
{
    double tini;
    int n_cell[NUM_DIM];
};


class c_simu
{
public:
    c_simu();
    ~c_simu();
    
    bool      stop;   // stop computation
    double    t, dt, tmax;    // time, timestep, max simulation time
    long int  it, itmax;      // time increment, max time increment
    c_grid    grid;           // simulation grid
    s_par     par;

    FILE*     fout;

    int initialise();   // load configurations
    int run();          // Runs the hydro evolution of the grid 
    // void prepareEvolution(long int it, double t, double dt); // computes wavespeeds,
    //     // updates boundaries
    // void evolve(long int it, double t, double dt);  // Time-integration
    // void setup();       // sets up the various hydro parameters before running
    // void updateCFL(long int it); // updates the value of the CFL depending on time-step
    // void timeScale();   // Computes the value of dt according to the CFL condition
    // void assignGrid(c_grid grid);
    // void printInfo();
    // void setDebugProfile();
    // void specialDebugProfile();

    // void setEvolvingRegion();   // discards uniform non-moving regions from evolution
    //                             // calculation

    // // in time.cpp (time integration methods)
    // #if TIME_INTEGRATION_ == RK1_
    //     void RK1Integration(long int it, double t, double dt);
    // #elif TIME_INTEGRATION_ == RK2_
    //     void RK2Integration(long int it, double t, double dt);
    // #elif TIME_INTEGRATION_ == RK3_
    //     void RK3Integration(long int it, double t, double dt);
    // #endif

    // IO methods
    // void dataDump(int it); 
};

#endif