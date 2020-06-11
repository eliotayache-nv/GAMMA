#ifndef FLUID_CELL_INTERFACE_H_
#define FLUID_CELL_INTERFACE_H_

//////////////////////////////////////////////////////////////////////////////////////////
#include "environment.h"
#include "fluid.h"

class c_cell;     // forward declaration

class c_cell_interface
{
public:
  c_cell_interface();
  ~c_cell_interface();

  // MEMBERS
  int  status; 
  int  memNumber;

  double x;               // position (in a single direction)
  double v;               // velocity (lab frame)
  double lfac;            // Lorentz factor (lab frame)
  double dim;             // orientation (orthogonal vector direction)
  double dA;              // surface area

  c_cell  *pCL;     // pointer to left cell
  c_cell  *pCR;     // pointer to right cells

  double F[NUM_Q];        // flux accross interface
  c_fluid_state  S;       // interface state
  c_fluid_state  SL;      // left state (after reconstruction)
  c_fluid_state  SR;      // right state (after reconstruction)

  void computeInterCellFlux(); // computes the flux vector accross the interface
};

#endif
