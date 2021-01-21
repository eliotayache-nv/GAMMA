/*
* @Author: eliotayache
* @Date:   2020-05-05 10:57:26
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-01-21 10:01:53
*/

#include "../environment.h"
#include "../grid.h"


void Grid::evolve(int it, double t, double dt){
  // only works in 2D for now

  copyState0();

  // first intermediate step
  update(dt);

  // second intermediate step
  updateGhosts(it, t);
  prepForUpdate(it, t);
  update(dt);
  // do not evolve border cells because they are going to be copied anyways
  // and it can lead to non-physical states
  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]+1; i <= iRbnd[j]-1; ++i){
      Cell *c = &Ctot[j][i];
      double dV  = c->G.dV;
      double dV0 = c->G0.dV;
      for (int q = 0; q < NUM_Q; ++q){
        double Q   = c->S.cons[q];
        double Q0  = c->S0.cons[q];
        c->S.cons[q] = 3./4. * Q0 * dV0 + 1./4. * Q * dV;
      }
    }
  }

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]; i <= iRbnd[j]-1; ++i){
      Interface *I = &Itot[j][i];
      double x  = I->x[MV];
      double x0 = I->x0[MV];
      I->x[MV] = 3./4. * x0 + 1./4. * x;
      I->computedA();
    }
  }

  CellGeomFromInterfacePos();

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]+1; i <= iRbnd[j]-1; ++i){
      Cell *c = &Ctot[j][i];
      for (int q = 0; q < NUM_Q; ++q){
        c->S.cons[q] /= c->G.dV;
      }
      c->S.cons2prim(c->G.x[r_]);
    }
  }
  
  // third intermediate step
  updateGhosts(it, t);
  prepForUpdate(it, t);
  update(dt);

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]+1; i <= iRbnd[j]-1; ++i){
      Cell *c = &Ctot[j][i];
      double dV  = c->G.dV;
      double dV0 = c->G0.dV;
      for (int q = 0; q < NUM_Q; ++q){
        double Q   = c->S.cons[q];
        double Q0  = c->S0.cons[q];
        c->S.cons[q] = 1./3. * Q0 * dV0 + 2./3. * Q * dV;
      }
    }
  }

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]; i <= iRbnd[j]-1; ++i){
      Interface *I = &Itot[j][i];
      double x  = I->x[MV];
      double x0 = I->x0[MV];
      I->x[MV] = 1./3. * x0 + 2./3. * x;
      I->computedA();
    }
  }

  CellGeomFromInterfacePos();

  #pragma omp parallel for default(shared)
  for (int j = jLbnd+1; j <= jRbnd-1; ++j){
    for (int i = iLbnd[j]+1; i <= iRbnd[j]-1; ++i){
      Cell *c = &Ctot[j][i];
      for (int q = 0; q < NUM_Q; ++q){
        c->S.cons[q] /= c->G.dV;
      }
      c->S.cons2prim(c->G.x[r_]);
    }
  }

}





























