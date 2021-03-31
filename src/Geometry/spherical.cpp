/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-03-31 11:20:39
*/


#include "../cell.h"
#include "../constants.h"

void Cell::computedV(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  double r = geom->x[r_];
  double th = geom->x[t_];
  double dr = geom->dx[r_];
  double dth = geom->dx[t_];
  double r1 = r - dr/2.;
  double r2 = r + dr/2.;
  // double th1 = th - dth/2.;
  // double th2 = th + dth/2.;
  double rsq = r1*r1 + r1*r2 + r2*r2;

  // geom->dV = fabs(2./3.*PI*(r2*r2*r2 - r1*r1*r1) * (cos(th1) - cos(th2)));
  // geom->dV = 2.*PI*r*r*fabs(sin(th))*dr*dth;

  // integrated but avoiding machine prec issues:
  geom->dV = 4./3.*PI*rsq*dr*fabs(sin(.5*dth)*sin(th));
    // fabs for negative thetas (should only exist for ghost cells, otherwise not allowed)

}


void Cell::computedl(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  double r   = geom->x[r_];
  double dr  = geom->dx[r_];
  double dth = geom->dx[t_];
  geom->dl[r_] = dr;
  geom->dl[t_] = r*dth;
  if (NUM_D == 3) exit(40);

}


void Cell::computeCentroid(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }
  
  double r = geom->x[r_];
  double dr = geom->dx[r_];
  double th = geom->x[t_];
  double dth = geom->dx[t_];
  double r2 = r*r;
  double dr2 = dr*dr;
  double th1 = th-dth/2.;
  double th2 = th+dth/2.;

  for (int d = 0; d < NUM_D; ++d){ geom->cen[d] = geom->x[d]; }
  // geom->cen[r_] = r + (2.* r * pow(dr, 2)) / (12. * pow(r, 2) + pow(dr, 2)); 
  geom->cen[r_] = 3.*r*(r2 + dr2) / (3*r2 + dr2); 
  geom->cen[t_] = (th1*cos(th1) - th2*cos(th2) + fabs(sin(th2)) - fabs(sin(th1))) 
    / (cos(th1) - cos(th2)); 

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
    // dA = 2.*PI*r*r*(cos(th-dth/2.) - cos(th+dth/2.));
    dA = 4.*PI*r*r*fabs(sin(.5*dth)*sin(th));
    // dA = 2*PI*r*r*fabs(sin(th))*dth;
  }
  else if (dim == t_){
    double dr = dx[0];
    double r1 = r - dr/2.;
    double r2 = r + dr/2.;
    double rsq = (r1+r2)*dr;
    // dA = PI*sin(th)*(r2*r2 - r1*r1);
    dA = PI*fabs(sin(th))*rsq;
      // fabs for negative thetas (should only exist for ghost cells, otherwise not allowed)
    // dA = 2.*PI*r*fabs(sin(th))*dr;
  }
  else {
    printf("This geometry is not implemented yet\n");
    exit(12);
  }

}


double Cell::regridVal(){

  double r = G.x[r_];
  double dr = G.dx[r_];
  double dth = G.dx[t_];
  double res = dr / (r*dth);
  user_regridVal(&res);
  return(res);

}

