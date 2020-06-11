#ifndef FLUID_CELL_H_
#define FLUID_CELL_H_

//////////////////////////////////////////////////////////////////////////////////////////
// STATUS FOR CELLS
enum{INACTIVE_STATUS_,ACTIVE_STATUS_,GHOST_STATUS_}

//////////////////////////////////////////////////////////////////////////////////////////
#include "environment.h"
#include "fluid.h"
#include "fluid_cell_interface.h"

struct s_cell_geometry
{
  double dV;                // volume
  double pos[NUM_DIM];      // position 
  double v[NUM_DIM];        // velocity (lab frame)
  double dl[NUM_DIM];       // width in every dimension
  double cen[NUM_DIM];      // centroid 
};

class c_fluid_cell
{
public:
  c_fluid_cell();
  c_fluid_cell(double rho, double v, double p);
  ~c_fluid_cell();

  int  status;
  int  memNumber;
  c_fluid_state   S;
  s_cell_geometry G;

  void update(double dt);
};

#endif













