/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:20
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 13:53:41
*/

#include "../fluid_cell.h"

void c_cell::computeCentroid()
{
    G.cen[x_] = G.pos[x_];
}