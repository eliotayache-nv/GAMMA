/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-21 07:10:10
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

    // printf("x y z\n");
    // for (int j = 0; j < grid.ncell[F1]; ++j)
    // {
    //   for (int i = 0; i < grid.nact[j]; ++i)
    //   {
    //     printf("%le %le %d\n", grid.C[j][i].G.x[x_], grid.C[j][i].G.x[y_], grid.C[j][i].nde_id);
    //   }
    // }

    dt = grid.prepForUpdate();
    // printf("%lf\n", dt);
    grid.update(dt);

    t += dt;
    it++;

    if (it == 100){ stop = true; }
  }

  // grid.print(PPP);
  // printf("x y z\n");
  // for (int j = 0; j < grid.ncell[F1]; ++j)
  // {
  //   for (int i = 0; i < grid.nact[j]; ++i)
  //   {
  //     printf("%le %le %le\n", grid.C[j][i].G.x[x_], grid.C[j][i].G.x[y_], grid.C[j][i].S.prim[RHO]);
  //   }
  // }

  return 0;

}
