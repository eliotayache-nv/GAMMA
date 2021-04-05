/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-03-31 11:25:30
*/


#include "../cell.h"
#include "../constants.h"

void Cell::computedV(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  double r = geom->x[r_];
  double dr = geom->dx[r_];
  geom->dV = 4.*PI*r*r*dr;

}


void Cell::computedl(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  double r   = geom->x[r_];
  double dr  = geom->dx[r_];
  geom->dl[r_] = dr;
  if (NUM_D == 2) exit(40);
  if (NUM_D == 3) exit(40);

}


void Cell::computeCentroid(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }\

  double r = geom->x[r_];
  double dr = geom->dx[r_];
  double r2 = r*r;
  double dr2 = dr*dr;

  for (int d = 0; d < NUM_D; ++d){ geom->cen[d] = geom->x[d]; }
  geom->cen[r_] = 3.*r*(r2 + dr2) / (3*r2 + dr2); 

}


void Cell::move(double xL, double xR){

  G.x[MV]  = (xR + xL) / 2.;
  G.dx[MV] = (xR - xL);
  computeAllGeom();

}


void Interface::computedA(){

  double r = x[r_];

  if (NUM_D == 3 or NUM_D == 2){
    exit(40);
  }

  if (dim == r_){
    dA = 4.*PI*r*r;
  }
  else {
    printf("This geometry is not implemented yet\n");
    exit(12);
  }

}


double Cell::regridVal(){

  double r = G.x[r_];
  double dr = G.dx[r_];
  double res = dr / r;
  user_regridVal(&res);
  return(res);

}

