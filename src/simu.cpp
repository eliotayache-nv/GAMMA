/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 13:38:45
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-16 14:58:29
*/
#include "simu.h"
#include "environment.h"

c_simu :: c_simu(){

  stop = false;
  t = 0;
  it = 0;
}

c_simu::~c_simu(){}

int c_simu :: initialise(){

  par.tini = 0.;
  t = par.tini;
  grid.initialise(par);
  return 0;
}

int c_simu :: run(){

    grid.print(RHO);
    return 0;
}
