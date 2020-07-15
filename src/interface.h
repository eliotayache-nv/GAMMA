#ifndef FLUID_CELL_INTERFACE_H_
#define FLUID_CELL_INTERFACE_H_

//////////////////////////////////////////////////////////////////////////////////////////
#include "environment.h"
#include "fluid.h"

class Cell;     // forward declaration

class Interface
{
public:
  Interface(int d=MV);
  ~Interface();

  // MEMBERS
  int  status; 
  int  memNumber;
  int  dim;               // orientation (orthogonal vector direction)

  double x[NUM_D];        // position (in a single direction)
  double v;               // velocity (only in MV dimension) (lab frame)
  double lfac;            // Lorentz factor (only in MV dimension) (lab frame)
  double dl[NUM_D-1];     // spatial extent
  double dA;              // surface area

  FluidState S,SL,SR;     
  double lL,lS,lR;        // wavespeeds
  double flux[NUM_Q];

  void wavespeedEstimates();
  void computeLambda();
  void computeFlux();
  FluidState starState(FluidState Sin, double lbda);

};

#endif
