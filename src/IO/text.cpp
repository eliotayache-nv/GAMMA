/*
* @Author: Eliot Ayache
* @Date:   2020-09-28 16:57:12
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-28 17:01:43
*/

#include "../simu.h"
#include "../mpisetup.h"
#include "../environment.h"
#include "../constants.h"

void Simu::reinitialise(){

  loadParams(&par);
  t = par.tini;
  grid.initialise(par);
  grid.initialGeometry();
  grid.initialValues();
  mpi_distribute(&grid);
  grid.prepForRun();

}