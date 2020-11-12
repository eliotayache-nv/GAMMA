/*
* @Author: Eliot Ayache
* @Date:   2020-06-12 11:54:57
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-11-09 15:45:54
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

  // printf("%le %le\n", v, dt);
  x[MV] += dt * v;
  computedA();

}
