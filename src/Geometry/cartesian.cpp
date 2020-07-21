/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-16 11:09:51
*/


#include "../cell.h"

void Cell::computedV(){

  G.dV = 1.;
  for (int d = 0; d < NUM_D; ++d){ G.dV *= G.dl[d]; }

}


void Cell::computeCentroid(){

  for (int d = 0; d < NUM_D; ++d){ G.cen[d] = G.x[d]; }

}


void Cell::move(double xL, double xR){

  G.x[MV]  = (xR + xL) / 2.;
  G.dl[MV] = (xR - xL);
  computedV();
  computeCentroid();

}