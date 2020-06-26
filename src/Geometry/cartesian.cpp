/*
* @Author: eliotayache
* @Date:   2020-05-05 15:17:31
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-25 08:51:36
*/


#include "../cell.h"

void Cell::computedV(){

    G.dV = 1.;
    for (int d = 0; d < NUM_D; ++d) G.dV *= G.dl[d];

}

void Cell::computeCentroid(){

    for (int d = 0; d < NUM_D; ++d) {
        G.cen[d] = G.x[d];
    }

}
