/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-04-11 20:07:40
*/


#include "../cell.h"

void Cell::computedV(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  geom->dV = 1.;
  for (int d = 0; d < NUM_D; ++d){ geom->dV *= geom->dx[d]; }

}

void Cell::computedl(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  for (int d = 0; d < NUM_D; ++d) geom->dl[d] = geom->dx[d];

}


void Cell::computeCentroid(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  for (int d = 0; d < NUM_D; ++d){ geom->cen[d] = geom->x[d]; }

}


void Cell::move(double xL, double xR){

  G.x[MV]  = (xR + xL) / 2.;
  G.dx[MV] = (xR - xL);
  computeAllGeom();

}


void Interface::computedA(){

  #if NUM_D == 1
    dA = 1;
  #else
    dA = dx[0];
  #endif

}


double Cell::regridVal(){

  double res = G.dx[MV];
  user_regridVal(&res);
  return(res);

}
