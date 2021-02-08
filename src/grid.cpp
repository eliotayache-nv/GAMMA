/*
* @Author: eliotayache
* @Date:   2020-05-06 09:26:35
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-02-08 21:35:23
*/

#include "grid.h"

Grid::Grid()
{
}

Grid::~Grid()
{

}


void Grid::prepForRun(){

  v2u();
  prim2cons();
  interfaceGeomFromCellPos();

}

void Grid::prepForUpdate(int it, double t){

  prim2cons();
  state2flux();
  computeNeighbors();
  movDir_ComputeLambda();   // has to be done before kinematics update (and flux calc)
  updateKinematics(it, t);
  computeFluxes();
  
}

