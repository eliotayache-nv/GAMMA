/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-15 16:26:21
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

    dt = grid.prepForUpdate();
    // printf("dt = %le\n", dt);
    grid.update(dt);

    t += dt;
    it++;

    if (it == 300){ stop = true; }
  }

  grid.print(RHO);
  // for (int j = 0; j < grid.nde_nax[F1]; ++j)
  // {
  //   for (int i = 0; i < grid.ntrack[j]; ++i)
  //   {
  //     printf("%le ", grid.Ctot[j][i].G.dV);
  //   }
  //   printf("\n");
  // }

  return 0;

}
