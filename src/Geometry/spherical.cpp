/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:20
* @Last Modified by:   eliotayache
* @Last Modified time: 2020-05-05 15:18:53
*/

#include "../fluid_cell.h"

void c_fluid_cell::computeCentroid()
{
    G.cen[x_] = G.pos[x_];
}