/*
* @Author: Eliot Ayache
* @Date:   2020-12-17 22:39:37
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-01-28 10:44:11
*/

#include "constants.h"

double Nalpha_, Nmp_, Nme_, Nqe_, NsigmaT_;

void normalizeConstants(double rhoNorm, double vNorm, double lNorm){

  Nalpha_ = alpha_ / (1. / (rhoNorm * vNorm * lNorm));
  Nmp_ = mp_ / (rhoNorm * pow(lNorm, 3));
  Nme_ = me_ / (rhoNorm * pow(lNorm, 3));
  Nqe_ = qe_ / (sqrt(rhoNorm) * pow(lNorm, 2.) * vNorm);
  NsigmaT_ =  sigmaT_ / pow(lNorm, 2.);

}
