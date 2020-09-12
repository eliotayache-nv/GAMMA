/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-12 12:48:21
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

    // auto start = std::chrono::high_resolution_clock::now();
    // grid.printCols(it);

    dt = grid.prepForUpdate(it, t);
    grid.update(dt);


    t += dt;
    it++;

    // auto finish = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = finish - start;
    // std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    if (it%500 == 0){ grid.printCols(it); }
    // if (it == 500){ stop = true; }
    if (t/c_ > 1.e10){ stop = true; }
  }

  return 0;

}
