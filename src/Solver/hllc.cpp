// /*
// * @Author: eliotayache
// * @Date:   2020-06-10 11:07:58
// * @Last Modified by:   eliotayache
// * @Last Modified time: 2020-11-11 09:51:07
// */

#include "../environment.h"
#include "../interface.h"


void Interface::computeFlux(){

  if (v < lL){
    S = SL;
    for (int q = 0; q < NUM_Q; ++q){ flux[q] = S.flux[dim][q] - v*S.cons[q]; }
                           // compensate for interface motion --^
  }                     
  if (lL <= v && v < lS){
    S = starState(SL, lL);
    for (int q = 0; q < NUM_Q; ++q){ 
      flux[q] = SL.flux[dim][q] + lL*(S.cons[q] - SL.cons[q]) - v*S.cons[q]; 
    }
  }
  if (lS <= v && v < lR){
    S = starState(SR, lR);
    for (int q = 0; q < NUM_Q; ++q){ 
      flux[q] = SR.flux[dim][q] + lR*(S.cons[q] - SR.cons[q]) - v*S.cons[q];
    }
  }
  if (lR <= v ){
    S = SR;
    for (int q = 0; q < NUM_Q; ++q){ flux[q] = S.flux[dim][q] - v*S.cons[q]; }
  }

  for (int q = 0; q < NUM_Q; ++q){ flux[q] *= dA; }


}
