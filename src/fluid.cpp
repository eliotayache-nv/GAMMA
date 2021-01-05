/*
* @Author: Eliot Ayache
* @Date:   2020-07-01 16:36:05
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-12-17 23:23:35
*/

#include "fluid.h"

double FluidState::lfac(){

  double u   = 0;
  double uu[NUM_D];

  for (int i = 0; i < NUM_D; ++i){
    uu[i] = prim[UU1+i];
    u += uu[i]*uu[i];
  }
  u = sqrt(u);

  return sqrt(1+u*u);

}

double FluidState::cs(){

  double u   = 0;
  double uu[NUM_D];

  for (int i = 0; i < NUM_D; ++i){
    uu[i] = prim[UU1+i];
    u += uu[i]*uu[i];
  }
  u = sqrt(u);

  return sqrt(1+u*u);

}