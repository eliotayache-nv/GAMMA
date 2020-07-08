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

  void prim2cons();
  void cons2prim(double pin = 0);
  void state2flux();

  double lfac();

};

#endif