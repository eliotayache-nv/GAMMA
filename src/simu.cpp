/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-30 16:46:58
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

  printf("Hello you! Welcome to GAMMA2D!\n");

  while (!stop){
    dt = grid.prepForUpdate(it, t);

    // printing grid (everything is ready right after grid prepare)
    if (it%50 == 0){ grid.printCols(it, t); }

    grid.update(dt);

    t += dt;
    it++;

    if ((worldrank == 0) and (it%1000 == 0)){ printf("it: %ld time: %le\n", it, t);}
    // if (it == 500){ stop = true; }
    if (t > 1.e10){ stop = true; }
  }

}
