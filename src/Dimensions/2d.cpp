/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-12 11:53:08
*/

#include "../environment.h"
#include "../grid.h"
#include "../interface.h"
#include "../cell.h"
#include "../array_tools.h"

static int splitGrid(int n_cell, int n_gst)
{
  /*
  int n_cell: total number of cells (non including ghost boundaries)
  int n_gst:  size of ghost regions.
  */
  int base = n_cell/worldsize;
  int rest = n_cell%worldsize;
  int nde_n_ax = base;
  if (worldrank<rest) nde_n_ax++;
  nde_n_ax+=2*n_gst;
  return nde_n_ax;
}

void c_grid :: initialise(s_par par)
{
  for (int d = 0; d < NUM_DIM; ++d) n_ax[d] = par.n_cell[d]+2*n_gst;
  nde_n_ax[MV_D_] = n_ax[MV_D_];
  nde_n_ax[FX_D1] = splitGrid(par.n_cell[FX_D1], n_gst);
  
  Ctot = array_2d<c_cell>(nde_n_ax[FX_D1],nde_n_ax[MV_D_]);
  I    = array_2d<c_interface>(nde_n_ax[FX_D1],nde_n_ax[MV_D_]+1);
  C    = &Ctot[n_gst];
}

void c_grid :: destruct()
{
  delete_array_2d(C);
  delete_array_2d(I);
}
