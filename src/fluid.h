#ifndef fluid_H_
#define fluid_H_ 

#include "environment.h"

class c_fluid_state
{
public:
  c_fluid_state(){}
  ~c_fluid_state(){}

  double prim[NUM_Q];
  double cons[NUM_Q];
  double flux[NUM_DIM][NUM_Q];

  void prim2cons();
  void cons2prim(double pin = 0);
  void state2flux();
};

#endif