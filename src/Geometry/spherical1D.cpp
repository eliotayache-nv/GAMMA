/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-02-02 15:31:52
*/


#include "../cell.h"
#include "../constants.h"

void Cell::computedV(){

  double r = G.x[r_];
  double dr = G.dx[r_];
  G.dV = 4.*PI*r*r*dr;

}


void Cell::computedl(){

  double r   = G.x[r_];
  double dr  = G.dx[r_];
  G.dl[r_] = dr;
  if (NUM_D == 2) exit(40);
  if (NUM_D == 3) exit(40);

}


void Cell::computeCentroid(){

  double r = G.x[r_];
  double dr = G.dx[r_];
  double r2 = r*r;
  double dr2 = dr*dr;

  for (int d = 0; d < NUM_D; ++d){ G.cen[d] = G.x[d]; }
  G.cen[r_] = 3.*r*(r2 + dr2) / (3*r2 + dr2); 

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

