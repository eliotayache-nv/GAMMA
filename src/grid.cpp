/*
* @Author: eliotayache
* @Date:   2020-05-06 09:26:35
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 18:19:44
*/

#include "grid.h"

c_grid :: c_grid()
{
}

c_grid::~c_grid()
{

}

void c_grid :: initialise()
{
    #if   NUM_DIM == 1
        C = array_1d<c_cell>(n_max)
        I = array_1d<c_cell>(n_max)
    #elif NUM_DIM == 2

    #elif NUM_DIM == 3

    #endif
}

void c_grid :: destruct()
{
    #if ndim_ == 1
        delete_array_1d(C);
        delete_array_1d(I);
    #elif ndim == 2
        delete_array_2d(C);
        delete_array_2d(I);
    #endif
}

