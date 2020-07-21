#include "cell.h"
#include "environment.h"

Cell::Cell(){

  for (int d = 0; d < NUM_D; ++d){
    for (int c = 0; c < NUM_Q; ++c){
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

  if (dim != MV){
    if (IR_lL != IL_lR) { dt_cand = a / (IR_lL - IL_lR); }
    else { dt_cand = 1.e15; }
  } 
  else {
    dt_cand = a / fmax(fabs(IL_lR), fabs(IR_lL));
  }

  if (dt_cand < 0) dt_loc = fmin(dt_loc, 1.e15);
  else dt_loc = fmin(dt_loc, dt_cand);

}


void Cell::update(double dt, double xL, double xR){


  for (int q = 0; q < NUM_Q; ++q){
    S.cons[q] *= G.dV;
  }

  for (int d = 0; d < NUM_D; ++d){
    for (int q = 0; q < NUM_Q; ++q){
      S.cons[q] += (flux[0][d][q] - flux[1][d][q]) * dt;
    }
  }

  move(xL, xR);

  for (int q = 0; q < NUM_Q; ++q){
    S.cons[q] /= G.dV;
  }

  S.cons2prim();

}












