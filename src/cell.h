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
  int    nde_id;          // id on node of cell (access with Ctot[0][id])
  int    nde_ind[NUM_D];  // location on grid (node-specific)
  double dt_loc = 1.e15;  // local max dt
  FluidState      S;
  s_cell_geometry G;
  double flux[2][NUM_D][NUM_Q];   // L and R fluxes in all dimensions and directions
  vector<int> neigh[NUM_D][2];    // neighboring cells in each direction (id)
  // double grad[NUM_D][NUM_Q];      // gradients in all dimensions

  void update(double dt, double xL, double xR);
  void move(double xL, double xR);
  void computedV();
  void computeCentroid();
  void computeAllGeom();

  void update_dt(int dim, double IL_lR, double IR_lL=0);

};

#endif













