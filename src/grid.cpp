/*
* @Author: eliotayache
* @Date:   2020-05-06 09:26:35
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-29 15:54:52
*/

#include "grid.h"

Grid::Grid()
{
}

Grid::~Grid()
{

}


void Grid::prepForRun(){

  // assignStatus();
  interfaceGeomFromCellPos();
  //TBC mode to come (ghost cells are not updated here yet)

}

void Grid::prepForUpdate(){

  updateGhosts();
  apply(&FluidState::prim2cons);
  apply(&FluidState::state2flux);

}

void Grid::update(){

  computeFluxes();

}
