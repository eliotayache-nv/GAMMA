/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   eliotayache
* @Last Modified time: 2020-11-11 09:51:28
*/
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

    // printing grid (everything is ready right after grid prepare)
    //if (it%20000 == 0){ grid.printCols(it, t); }
    //if (it%10000 == 0){ grid.printCols(it, t); }
    // if (it%2000 == 0){ grid.printCols(it, t); }
    if (it%1000 == 0){ grid.printCols(it, t); }

    //if ((worldrank == 0) and (it%1000 == 0)){ printf("it: %ld time: %le\n", it, t);}
    if ((worldrank == 0) and (it%100 == 0)){ printf("it: %ld time: %le\n", it, t);}

    //if (it > 8870000){ stop = true; }
    //if (it > 50000){ stop = true; }
    //if (t > 2.5e6){ stop = true; }
    // if (it > 1000){ stop = true; }
    if (t > 1.e7){ stop = true; }
  }

}
