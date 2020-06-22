/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:20
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-18 11:28:16
*/

#include "../fluid_cell.h"

void Cell::computeCentroid()
{
    G.cen[x_] = G.x[x_];
}