#include "cell.h"
#include "environment.h"

Cell::Cell(){

  for (int d = 0; d < NUM_D; ++d){
    for (int c = 0; c < NUM_Q; ++c){
      flux[0][d][c] = 0;
      flux[1][d][c] = 0;
    }
  }

  for (int d = 0; d < NUM_D; ++d){ S.prim[VV1+d] = 0.; }

  for (int q = 0; q < NUM_Q; ++q){
    grad_mv[q] = 0;
    grad_f1[q] = 0;
  }

  #if SHOCK_DETECTION_ == ENABLED_
    isShocked = false;
    Sd = -1.;
    pspec = 0;
  #endif

}

Cell::~Cell(){}

void Cell::computeAllGeom(s_cell_geometry *geom){

  if (geom==NULL) {geom = &G; }

  computedl(geom);
  computedV(geom);
  computeCentroid(geom);

}

void Cell::resetLocaldt(){

  dt_loc = 1.e15;

}


void Cell::update_dt(int dim, Interface IL, Interface IR){

  double a = G.dl[dim];
  double dt_cand;

  if (dim == MV){

    double l, v;
    double dt_candL, dt_candR;
    // wave from left side:
    l = IL.lR;
    v = IR.v;
    if (l > v) { dt_candL = a / (l - v); }
    else { dt_candL = 1.e15; }

    // wave from right side:
    l = IR.lL;
    v = IL.v;
    if (l < v) { dt_candR = a / (v - l); }
    else { dt_candR = 1.e15; }

    dt_cand = fmin(dt_candL, dt_candR);

    if (dt_cand < 0){
      printf("MV %le %le %le %le\n", dt_cand, a, l, v);
      exit(30);
    }
  } 

  else {

    double l;
    double xI = IL.x[dim];
    double xC = G.x[dim];

    if (xI < xC){
      l = IL.lR;
      if (l > 0) { dt_cand = a/l;  } else { dt_cand = 1.e15; }
    }
    else{
      l = -IL.lL;
      if (l > 0) { dt_cand = a/l;  } else { dt_cand = 1.e15; }      
    }
    if (dt_cand < 0){
      printf("F1 %le %le %le\n", dt_cand, a, l);
      exit(31);
    }
  }


  dt_loc = fmin(dt_loc, dt_cand);

}


void Cell::update(double dt, double xL, double xR){


  if (isnan(flux[1][1][TR1])){
    printf("%d %d\n", nde_ind[F1], nde_ind[MV]);
  }

  for (int q = 0; q < NUM_Q; ++q){
    S.cons[q] *= G.dV;
  }

  for (int d = 0; d < NUM_D; ++d){
    for (int q = 0; q < NUM_Q; ++q){
      S.cons[q] += (flux[0][d][q] - flux[1][d][q]) * dt;
    }
  }
  sourceTerms(dt); userSourceTerms(dt); 

  #if LOCAL_SYNCHROTRON_ == ENABLED_
    radiativeSourceTerms(dt);
  #endif
    
  move(xL, xR);

  for (int q = 0; q < NUM_Q; ++q){
    S.cons[q] /= G.dV;
  }
  S.cons2prim(G.x[r_]);

  if (isnan(S.prim[TR1])){
    printf("%le\n", flux[1][1][TR1]);
    exit(30);
  }

}












