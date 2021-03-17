/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-03-17 19:25:33
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

  #if LOCAL_SYNCHROTRON_ == ENABLED_
    grid.apply(&Cell::radiation_initGammas);
  #endif

}


void Simu::run(){
  
  grid.updateGhosts(it, t);
  grid.printCols(it, t);

  while (!stop){

    grid.regrid();

    #if SHOCK_DETECTION_ == ENABLED_
      grid.apply(&Cell::resetShock);
    #endif
      
    grid.updateGhosts(it, t);
    grid.prepForUpdate(it, t);

    #if SHOCK_DETECTION_ == ENABLED_
      grid.apply(&Cell::detectShock);
    #endif

    dt = CFL_ * grid.collect_dt();

    grid.evolve(it, t, dt);

    t += dt;
    it++;

    #if LOCAL_SYNCHROTRON_ == ENABLED_
      grid.apply(&Cell::radiation_injectParticles);
    #endif

    dataDump();
    runInfo();
    evalEnd();
    
  }

}
