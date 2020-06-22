/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-18 15:58:22
*/
#include "simu.h"
#include "mpisetup.h"
#include "environment.h"

Simu :: Simu(){

  stop = false;
  t = 0;
  it = 0;

}

Simu::~Simu(){}

int Simu :: initialise(){

  int status = 0;
  loadParams(&par);
  t = par.tini;
  grid.initialise(par);
  status = grid.initialGeometry();
  status = grid.initialValues();
  mpi_distribute(&grid);

  return status;

}

int Simu :: run(){

    grid.print(RHO);
    grid.destruct();
    return 0;
}
