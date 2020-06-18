#ifndef FLUID_CELL_H_
#define FLUID_CELL_H_

//////////////////////////////////////////////////////////////////////////////////////////
// STATUS FOR CELLS
enum{INACTIVE_STATUS_,ACTIVE_STATUS_,GHOST_STATUS_};

//////////////////////////////////////////////////////////////////////////////////////////
#include "environment.h"
#include "fluid.h"
#include "interface.h"

struct s_cell_geometry{

  double dV;                // volume
  double pos[NUM_D];      // position 
  double v[NUM_D];        // velocity (lab frame)
  double dl[NUM_D];       // width in every dimension
  double cen[NUM_D];      // centroid 

};

class c_cell{

public:
  c_cell();
  ~c_cell();

  int  status;
  int  memNumber;
  c_fluid_state   S;
  s_cell_geometry G;

  void update(double dt);

};

#endif













