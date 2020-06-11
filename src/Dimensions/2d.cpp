/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:58:15
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 19:27:47
*/

#include "grid.h"
#include "environment.h"

void c_grid :: initialise(s_par par)
{
  for (int d = 0; d < NUM_DIM; ++d) n_ax[d] = par.n_cell[d]+2*n_gst;
  nde_n_ax[MOV_DIM] = n_ax[MOV_DIM];

  C = array_2d<c_cell>(n_max)
  I = array_2d<c_cell>(n_max)
}


void c_grid :: destruct()
{
  delete_array_2d(C);
  delete_array_2d(I);
}
