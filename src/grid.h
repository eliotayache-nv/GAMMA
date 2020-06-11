#ifndef FLUID_GRID_H_
#define FLUID_GRID_H_

#include "environment.h"
#include "cell.h"
#include "array_tools.h"


//////////////////////////////////////////////////////////////////////////////////////////
struct s_config_setup;

class c_grid
{
public:
  c_grid();
  ~c_grid();

  // total grid info
  int n_ax[NUM_DIM];  // number of cells in each direction to be allocated
                        // n_ax[MOV_DIM] = n_max 
  int n_max;          // max number of cells on moving dim
  int n_gst;          // number of ghost cells on each side of the grid 

  // node-specific info
  int nde_n_ax[NUM_DIM];  // number of cells in each direction

  #if   NUM_DIM == 1
    c_cell        *C;
    c_cell_interface    *I;     // moving interfaces

    int n_act;                  // number of active cells in track
    int iLactive, iRactive;     // indexes of leftmost and rightmost active cells on track

  #elif NUM_DIM == 2
    c_cell        **C;
    c_cell_interface    **I;

    int *n_act;
    int *iLactive, *iRactive;

  #elif NUM_DIM == 3
    c_cell        ***C;
    c_cell_interface    ***I;

    int **n_act;
    int **iLactive, **iRactive;

  #endif


  // bool initialised;

  // int ncells_active;          // number of active cells at a given moment (ACTIVE)
  // int ncells_total;           // total number of cells (ghost, active)
  // int ncells_max;             // max number of cells (ghost, active, inactive)
  // int iLactive, iRactive;     // indexes of leftmost and rightmost active cells
  // int iLlim, iRlim;           // indexes of left and right evolving region limits
  // c_cell        *C;
  // c_cell_interface    *I;

  // // METHODS
  void initialise();          // allocate memory for maximum number of cells
  void destruct();            // free memory
  // void loadInitialConfig(s_config_setup config);   // Fill grid with initial conditions
  // void reloadConfig(s_config_setup *pconfig); // Fill grid with result of previous run
  // void updateBoundaries();    // sets up the boundary conditions

  // // in geometry.h (geometry setup and update)
  // void setupGeometry(s_config_setup config);  // compute cell geometrical variables
  // void geomFromInterfacePositions();  // new geom for new interface positions
  //                   // analog to las bit of moveMesh()

  // // in hydro.cpp (hydro evolution)
  // void updateWaveSpeeds();    // Updates the wavespeeds at all interfaces
  // void updateCells(long int it, double t, double dt);// Computes the next state of the grid
  // void updateInterCellFluxes();// Computes the flux accross all interfaces of the grid
  // void reconstructInterfaceStates(double dt); // Runs reconstruction for all interfaces in grid
  // void updateSolver(long int it, double t, double dt); // chooses the solver to be used
  //                            // for each interface
  // void forceSolverToHLLC();   // Forces solver to HLLC on all interfaces. (Useful for
  //               // ALE_TYPE_ == CD_)

  // // in time.cpp (time integration methods)
  // void firstOrderUpdate(long int it, double t, double dt); // updates cells and moves mesh

  // // in hydro_moving_mesh.cpp (moving mesh methods)
  // void updateMeshKinematics(long int it, double t); // updates the Mesh kinematics in accordance
  //                     // with grid states
  // void moveMesh(double t, double dt); // Moves the mesh according to kinematics 
  //                   // specified before
  // void regrid(s_config_setup config, int boolRefine, double t);  // merges cells that are too small.
  //                    // splits cells that are too big.
  // void flagGradients(); // computes second derivative of velocity over all grid.
  // void refine(s_config_setup config, int ix);     // refines for index ix
  // void derefine(s_config_setup config, int ix);   // de-refines for index ix
  // void findMinInactive();   // updates minInactiveC and pFirstInactiveI
  // void correctCellsForMovingMesh(double dt);  // Corrects cell values to account for
  //                       // interfaces position updates.
  // void sort(); // sorts the cells in memory to be ordered again
  // void moveRightBoundary(double dt); // moves the right side boundary when activated

  // #if REFINEMENT_CRITERION_ == USER_REFINEMENT_
  // void userRefine(s_config_setup, int ix);        // custom refinement decision
  // #endif
  // #if DEREFINEMENT_CRITERION_ == USER_DEREFINEMENT_
  // void userDerefine(s_config_setup, int ix);      // custom derefinement decision
  // #endif
  

  // // extending fluid functions to the grid
  // void prim2cons();
  // void prim2aux();
  // void cons2prim();
  // void state2flux();
  // void cons2all(long int it, double t, double dt);


  // // RK3 methods
  // void RK3CopyState0();
  // void RK3ComputeState2(long int it, double t, double dt);
  // void RK3ComputeState3(long int it, double t, double dt);

};

#endif
