#ifndef fluid_H_
#define fluid_H_ 

#include "environment.h"

class FluidState{

public:
  FluidState();
  ~FluidState();

  double prim[NUM_Q];
  double cons[NUM_Q];
  double flux[NUM_D][NUM_Q];

  void prim2cons(double r);
  void cons2prim(double r, double pin = 0);
  void cons2prim_user(double *rho, double *p, double *u);
  void state2flux(double r);

  double lfac();
  double cs();

};

#endif