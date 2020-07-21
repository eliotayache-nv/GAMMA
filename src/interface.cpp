/*
* @Author: Eliot Ayache
* @Date:   2020-06-12 11:54:57
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-16 10:40:46
*/

#include "interface.h"

Interface::Interface(int d)
{
  dim = d;

  lL = 0.;
  lS = -1.;   // stays at -1 if unused.
  lR = 0.;
}

Interface::~Interface()
{
  
}

void Interface::move(double dt){

  x[MV] += dt * v;
  // TBC add dl update

}
