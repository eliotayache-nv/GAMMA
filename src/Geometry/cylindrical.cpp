/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-10-11 20:42:36
*/


#include "../cell.h"
#include "../constants.h"

void Cell::computedV(){

  double r = G.x[r_];
  double dr = G.dx[r_];
  double dz = G.dx[zcyl_];
  double r1 = r - dr/2.;
  double r2 = r + dr/2.;

  G.dV = PI*(r2*r2 - r1*r1) * dz;

}


void Cell::computedl(){

  double dr = G.dx[r_];
  double dz = G.dx[zcyl_];

  G.dl[r_] = dr;
  G.dl[zcyl_] = dz;

  if (NUM_D == 3) exit(40);

}


void Cell::computeCentroid(){

  double r = G.x[r_];
  double dr = G.dx[r_];
  double r1 = r - dr;
  double r2 = r + dr;

  for (int d = 0; d < NUM_D; ++d){ G.cen[d] = G.x[d]; }
  G.cen[r_] = 2./3. * (r2*r2*r2 - r1*r1*r1) / (r*dr);   // rel notes

}


void Cell::move(double xL, double xR){

  G.x[MV]  = (xR + xL) / 2.;
  G.dx[MV] = (xR - xL);
  computeAllGeom();

}


void Interface::computedA(){

  double r = x[r_];

  if (NUM_D == 3){
    printf("3D not implemented yet!\n");
    exit(12);
  }

  if (dim == r_){
    double dz = dx[0];
    dA = 2.*PI*r*dz;
  }
  else if (dim == zcyl_){
    double dr = dx[0];
    double r1 = r - dr/2.;
    double r2 = r + dr/2.;

    dA = PI*(r2*r2 - r1*r1);
  }
  else {
    printf("This geometry is not implemented yet\n");
    exit(12);
  }

}


double Cell::regridVal(){

  double res = G.dx[zcyl_] / G.dx[r_];
  user_regridVal(&res);
  return(res);

}

