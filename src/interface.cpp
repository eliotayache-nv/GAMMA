/*
* @Author: Eliot Ayache
* @Date:   2020-06-12 11:54:57
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-01 16:25:20
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
