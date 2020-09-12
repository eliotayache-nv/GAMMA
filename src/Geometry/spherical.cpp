/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-12 00:02:29
*/


#include "../cell.h"
#include "../constants.h"

void Cell::computedV(){

  double r = G.x[r_];
  double th = G.x[t_];
  double dr = G.dx[r_];
  double dth = G.dx[t_];
  double r1 = r - dr/2.;
  double r2 = r + dr/2.;
  double th1 = th - dth/2.;
  double th2 = th + dth/2.;

  G.dV = fabs(2./3.*PI*(r2*r2*r2 - r1*r1*r1) * (cos(th1) - cos(th2)));
  // G.dV = 2.*PI*r*r*fabs(sin(th))*dr*dth;
    // fabs for negative thetas (should only exist for ghost cells, otherwise not allowed)

}


void Cell::computedl(){

  double r   = G.x[r_];
  double dr  = G.dx[r_];
  double dth = G.dx[t_];
  G.dl[r_] = dr;
  G.dl[t_] = r*dth;
  if (NUM_D == 3) exit(40);

}


void Cell::computeCentroid(){

  double r = G.x[r_];
  double dr = G.dx[r_];

  G.cen[r_] = r + (2.* r * pow(dr, 2)) / (12. * pow(r, 2) + pow(dr, 2)); 
  // for (int d = 0; d < NUM_D; ++d){ G.cen[d] = G.x[d]; }

}


void Cell::move(double xL, double xR){

  G.x[MV]  = (xR + xL) / 2.;
  G.dx[MV] = (xR - xL);
  computeAllGeom();

}


void Interface::computedA(){

  double r = x[r_];
  double th = x[t_];

  if (NUM_D == 3){
    printf("3D not implemented yet!\n");
    exit(12);
  }

  if (dim == r_){
    double dth = dx[0];
    dA = 2.*PI*r*r*(cos(th-dth/2.) - cos(th+dth/2.));
    // dA = 2*PI*r*r*fabs(sin(th))*dth;
  }
  else if (dim == t_){
    double dr = dx[0];
    double r1 = r - dr/2.;
    double r2 = r + dr/2.;

    dA = PI*sin(th)*(r2*r2 - r1*r1);
      // fabs for negative thetas (should only exist for ghost cells, otherwise not allowed)
    // dA = 2.*PI*r*fabs(sin(th))*dr;
  }
  else {
    printf("This geometry is not implemented yet\n");
    exit(12);
  }

}


double Cell::regridVal(){

  return(G.x[r_] * G.dx[r_]);

}

