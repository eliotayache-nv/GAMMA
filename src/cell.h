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
  double dx[NUM_D];       // width in every dimension
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
  FluidState      S, S0;  // 0 for higher time-stepping
  s_cell_geometry G, G0;  // 0 for higher time-stepping
  double flux[2][NUM_D][NUM_Q];   // L and R fluxes in all dimensions and directions
  vector<int> neigh[NUM_D][2];    // neighboring cells in each direction (id)
  double grad_mv[NUM_Q];          // gradient in moving direction

  #if SHOCK_DETECTION_ == ENABLED_
    bool    isShocked;
    double  Sd;
    void detectShock();
    void resetShock();
  #endif

  #if LOCAL_SYNCHROTRON_ == ENABLED_
    void radiation_initGammas();
    void radiation_injectParticles();
    void radiation_apply_trac2gammae();
    void radiativeSourceTerms(double dt);
  #endif

  void update(double dt, double xL, double xR);
  void move(double xL, double xR);
  void sourceTerms(double dt);
  void userSourceTerms(double dt);
  void computedV();
  void computedl();
  void computeCentroid();
  void computeAllGeom();

  void resetLocaldt();
  void update_dt(int dim, Interface IL, Interface IR=NULL);

  // AMR
  double regridVal();
  void user_regridVal(double *res);

};

#endif













