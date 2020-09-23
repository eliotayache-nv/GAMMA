#ifndef FLUID_GRID_H_
#define FLUID_GRID_H_

#include "environment.h"
#include "cell.h"

//////////////////////////////////////////////////////////////////////////////////////////
class Grid{

  public:
  Grid();
  ~Grid();

  // total grid info
  int nsimu;         // number of cells in simu
  int nax[NUM_D];    // number of cells in each direction to be allocated
                        // nax[MV] = n_max (includes ghosts)
  int ncell[NUM_D];  // number of cells in simulation domain
  int ngst;          // number of ghost cells on each side of the grid

  // node-specific info
  int nde_ntot;        // total number of allocated cells on node
  int nde_nax[NUM_D];  // number of cells in each direction (including ghosts & inactive)
  int nde_ncell[NUM_D];  // initial number of cells in each direction 
  int origin[NUM_D];    // coordinates of C[0][0] in simulation domain

  #if   NUM_D == 1
    Cell         *C;
    Interface    *I;     // moving interfaces
    int n_act;           // adaptive number of active cells in track

  #elif NUM_D == 2
    Cell         **Cinit; // initial simulation domain
    Cell         **C;     // physical grid (without ghost tracks, but ghosts cells in MV)
    Cell         **Ctot;  // grid with ghost tracks
    Interface    **I;     // moving dim interfaces (without ghost to ghost in moving dim)
                            // does include ghost tracks
    Interface    **Itot;  // moving dim interfaces 
    int *nact;    // adaptive number of active cells in track
    int *ntrack;  // adaptive number of active and ghost cells in track
    int jLbnd, jRbnd, *iLbnd, *iRbnd; // tot index of first ghost cell out of active grid
    int *ibig, *ismall; // indexes of biggest and smallest cells along track

  #elif NUM_D == 3
    Cell         ***C;
    Interface    ***I;
    int **n_act;
    int **iLactive, **iRactive;

  #endif


  // METHODS
  // initialisation
  void initialise(s_par par);          // allocate memory for maximum number of cells
  int  initialGeometry();
  int  initialValues();
  void destruct();            // free memory
  void print(int var);
  void printCols(int it);
  void interfaceGeomFromCellPos();
  void interfaceGeomFromCellPos(int j); // only one track

  // update
  void prepForRun();
  void computeNeighbors(bool print=false);
  void assignId(int ind[NUM_D]);
  void movDir_ComputeLambda();
  void updateKinematics();
  void userKinematics();
  void gradients(Cell *c);
  double prepForUpdate(int it, double t);
  void computeFluxes();
  double collect_dt();
  void updateGhosts(int it, double t);
  void userBoundaries(int it, double t);
  void reconstructStates( int j, int i, int dim, 
                          int idn=-1,
                          Interface *Int = NULL );
  void update(double dt);

  // AMR
  void regrid();
  void targetRegridVictims( int j);
   int checkCellForRegrid(int j, int i);
  void applyRegrid(int j, int i, int action);
  void split(int j, int i);
  void merge(int j, int i);

  // toools
  void prim2cons();
  void state2flux();
  void apply(void (Cell::*func)());
  void apply(void (FluidState::*func)(), bool noborder=false);
    // overloading for FluidState

  void mpi_exchangeGhostTracks();

};

#endif
