#include "simu.h"
#include "mpisetup.h"
#include "environment.h"
#include "constants.h"

Simu::Simu(){

  stop = false;
  t = 0;
  it = 0;

}


Simu::~Simu(){}


void Simu::initialise(){

  loadParams(&par);
  t = par.tini;
  grid.initialise(par);
  grid.initialGeometry();
  grid.initialValues();
  mpi_distribute(&grid);
  grid.prepForRun();

}


void Simu::run(){
  
  grid.printCols(it, t);

  while (!stop){


    grid.regrid();

    grid.updateGhosts(it, t);
    grid.prepForUpdate(it, t);

    dt = CFL_ * grid.collect_dt();

    grid.evolve(it, t, dt);

    t += dt;
    it++;

    if (it%2000 == 0){ grid.printCols(it, t); }
    // if (it%100 == 0){ grid.printCols(it, t); }
    // if (it%10 == 0){ grid.printCols(it, t); }

    double t_C = 3.e3;
    double t_S = 9.e4;
    double t_g = 6.498e6;
    double t_M = 2.019e7;

    if ((t>=t_C) and ((t-dt)<t_C)){grid.printCols(it,t);}
    if ((t>=t_S) and ((t-dt)<t_S)){grid.printCols(it,t);}
    if ((t>=t_g) and ((t-dt)<t_g)){grid.printCols(it,t);}
    if ((t>=t_M) and ((t-dt)<t_M)){grid.printCols(it,t);}

    if ((worldrank == 0) and (it%100 == 0)){ printf("it: %ld time: %le\n", it, t);}

    if (it > 5000000){ stop = true; }
    // if (it > 100){ stop = true; }
    // if (t > 4.e3){ stop = true; }
    if (t > 5.e7){ stop = true; }
  }

}
