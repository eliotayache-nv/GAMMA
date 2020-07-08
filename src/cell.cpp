#include "cell.h"
#include "environment.h"

Cell::Cell(){

  for (int d = 0; d < NUM_D; ++d){
    for (int c = 0; c < NUM_C; ++c){
      flux[0][d][c] = 0;
      flux[1][d][c] = 0;
    }
  }

}

Cell::~Cell(){}

void Cell::computeAllGeom(){

  computedV();
  computeCentroid();

}

void Cell::update_dt(int dim, double IL_lR, double IR_lL){

  double a = G.dl[dim];
  double dt_cand;

  if (dim != MV)
  {
    if (IR_lL != IL_lR) { dt_cand = a / (IR_lL - IL_lR); }
    else { dt_cand = 1.e15; }
  } 
  else {
    dt_cand = a / fmax(IL_lR, IR_lL);
  }

  if (dt_cand < 0) dt_loc = fmin(dt_loc, 1.e15);
  else dt_loc = fmin(dt_loc, dt_cand);

}













