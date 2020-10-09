/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-09 12:57:58
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
    if (it%1 == 0){ grid.printCols(it, t); }

    if ((worldrank == 0) and (it%10 == 0)){ printf("it: %ld time: %le\n", it, t);}
    if (it > 1000){ stop = true; }
    if (t > 1.e10){ stop = true; }
  }

}
