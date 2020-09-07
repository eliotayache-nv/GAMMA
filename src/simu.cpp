/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-07 16:47:00
*/
#include "simu.h"
#include "mpisetup.h"
#include "environment.h"

Simu::Simu(){

  stop = false;
  t = 0;
  it = 0;

}

Simu::~Simu(){}

int Simu::initialise(){

  int status = 0;

  loadParams(&par);
  t = par.tini;
  grid.initialise(par);
  status = grid.initialGeometry();
  status = grid.initialValues();
  mpi_distribute(&grid);
  grid.prepForRun();

  return status;

}

int Simu::run(){

  while (!stop){
    dt = grid.prepForUpdate(it);
    grid.update(dt);

    t += dt;
    it++;
    if (it == 100){ stop = true; }
  }

  grid.printCols();

  return 0;

}
