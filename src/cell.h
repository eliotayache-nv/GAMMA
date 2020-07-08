#ifndef FLUID_CELL_H_
#define FLUID_CELL_H_

//////////////////////////////////////////////////////////////////////////////////////////
// STATUS FOR CELLS
enum{INACTIVE_STATUS_,ACTIVE_STATUS_,BOUND_STATUS_,GHOST_STATUS_};

//////////////////////////////////////////////////////////////////////////////////////////
#include "environment.h"
#include "fluid.h"
#include "interface.h"

struct s_cell_geometry{

  double dV;              // volume
  double x[NUM_D];        // position 
  double dl[NUM_D];       // width in every dimension
  double cen[NUM_D];      // centroid 

};

class Cell{

public:
  Cell();
  ~Cell();

  int    status;
  int    memNumber;
  double dt_loc = 1.e15;  // local max dt
  FluidState      S;
  s_cell_geometry G;
  double flux[2][NUM_D][NUM_C];  // L and R fluxes in all dimensions and directions

  void update(double dt);
  void computedV();
  void computeCentroid();
  void computeAllGeom();

  void update_dt(int dim, double IL_lR, double IR_lL=0);

};

#endif













