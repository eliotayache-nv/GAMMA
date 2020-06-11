/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 14:14:51
*/


#include "../fluid_cell.h"

void c_cell::computeCentroid()
{
    G.cen[x_] = G.pos[x_]
        + (2.* G.pos[x_] * pow(G.dl[x_], 2)) / (12. * pow(G.pos[x_], 2) + pow(G.dl[x_], 2)); 
}
