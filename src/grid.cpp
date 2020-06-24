/*
* @Author: eliotayache
* @Date:   2020-05-06 09:26:35
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-24 16:19:27
*/

#include "grid.h"

Grid :: Grid()
{
}

Grid::~Grid()
{

}

void Grid::prepEvol(){

  apply(&FluidState::prim2cons);
  apply(&FluidState::state2flux);

}