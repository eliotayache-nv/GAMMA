/*
* @Author: eliotayache
* @Date:   2020-05-06 09:26:35
* @Last Modified by:   eliotayache
* @Last Modified time: 2020-06-11 12:48:30
*/

#include "fluid_grid.h"

c_fluid_grid :: c_fluid_grid()
{
    initialized = false;
}

c_fluid_grid::~c_fluid_grid()
{

}

c_fluid_grid :: initialize()
{
    #if   NUM_DIM == 1
        // C = array_1d<c_fluid_cell>(n_max)
        // I = array_1d<c_fluid_cell>(n_max)
    #elif NUM_DIM == 2

    #elif NUM_DIM == 3

    #endif
}

void c_fluid_grid :: destruct()
{
    #if ndim_ == 1
        delete_array_1d(C);
        delete_array_1d(I);
    #elif ndim == 2
        delete_array_2d(C);
        delete_array_2d(I);
    #endif
    initialized = false;
}

