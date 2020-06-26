#ifndef FLUID_CELL_INTERFACE_H_
#define FLUID_CELL_INTERFACE_H_

//////////////////////////////////////////////////////////////////////////////////////////
#include "environment.h"
#include "fluid.h"

class Cell;     // forward declaration

class Interface
{
public:
  Interface();
  ~Interface();

  // MEMBERS
  int  status; 
  int  memNumber;

  double x[NUM_D];        // position (in a single direction)
  double v;               // velocity (only in MV dimension) (lab frame)
  double lfac;            // Lorentz factor (only in MV dimension) (lab frame)
  double dim;             // orientation (orthogonal vector direction)
  double dl[NUM_D-1];     // spatial extent
  double dA;              // surface area

  // Cell  *pCL;     // pointer to left cell
  // Cell  *pCR;     // pointer to right cells

  // double F[NUM_Q];        // flux accross interface
  // FluidState  S;       // interface state
  // FluidState  SL;      // left state (after reconstruction)
  // FluidState  SR;      // right state (after reconstruction)

  // void computeInterCellFlux(); // computes the flux vector accross the interface
};

#endif
