/*
* @Author: eliotayache
* @Date:   2020-05-05 10:31:06
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-04-11 19:59:25
*/

#include "../../environment.h"
#include "../../grid.h"
#include "../../simu.h"
#include "../../constants.h"

static double xmin = -0.35;
static double xmax = 1.;
static double ymin = -1.;
static double ymax = 1.;

static double Nx = 3160;
static double rho_ref = 1.;
static double v_ref   = 0.;
static double p_ref   = 1.e2;
static double alpha   = 1.;
static double L       = 0.3;
static double angle   = 0;//M_PI/4.;


void loadParams(s_par *par){

  par->tini      = 0.;
  par->ncell[x_] = Nx;
  par->nmax      = Nx+50;    // max number of cells in MV direction
  par->ngst      = 2;

}

int Grid::initialGeometry(){

  double x = xmax - xmin;
  double y = ymax - ymin;
  for (int i = 0; i < ncell[x_]; ++i){
    Cell *c = &Cinit[i];
    c->G.x[x_]  = (double) x*(i+0.5)/ncell[x_] + xmin;
    c->G.dx[x_] =          x/ncell[x_];

    c->computeAllGeom();
  }
  return 0;

}


static double Jp(double v, double cs)
{
  return(1./2. * log((1+v)/(1-v)) \
    + 1./sqrt(GAMMA_-1.) * log((sqrt(GAMMA_-1)+cs) / (sqrt(GAMMA_-1)-cs)));
}

static double Jm(double v, double cs)
{
  return(1./2. * log((1+v)/(1-v)) \
    - 1./sqrt(GAMMA_-1.) * log((sqrt(GAMMA_-1)+cs) / (sqrt(GAMMA_-1)-cs)));
}


int Grid::initialValues(){

  double h_ref  = 1. + p_ref * GAMMA_ / (GAMMA_ - 1.) / rho_ref; 
  double cs_ref = sqrt(GAMMA_ * p_ref / (rho_ref*h_ref));

  // Riemann invariants
  // double lp_ref = (v_ref + cs_ref) / (1 + v_ref*cs_ref);
  // double Jp_ref = Jp(v_ref, cs_ref);
  double Jm_ref = Jm(v_ref, cs_ref);

  for (int j = 0; j < ncell[F1]; ++j){
    for (int i = 0; i < ncell[MV]; ++i){
      Cell *c = &Cinit[i];
      double x = c->G.x[x_];
      double y = c->G.x[y_];
      double x_eff = x*cos(angle) - y*sin(angle); // rotation
      // double y_eff = x*sin(angle) + y*cos(angle);

      double rho, v, p;

      if (fabs(x_eff) < L and fabs(y < 0.3)){
        rho = rho_ref * (1. + alpha * pow((x_eff/L)*(x_eff/L) - 1. ,4.));
        p = p_ref * pow(rho / rho_ref, GAMMA_);
        double h = 1. + p * GAMMA_ / (GAMMA_ - 1.) / rho; 
        double cs = sqrt(GAMMA_ * p / (rho*h));        
        double temp = exp(2.*Jm_ref)
          * pow(((sqrt(GAMMA_-1.)+cs) / (sqrt(GAMMA_-1.)-cs)),(2./sqrt(GAMMA_-1.)));
        v = (temp - 1.)/(temp + 1.);

        c->S.prim[RHO] = rho;
        c->S.prim[PPP] = p;
        c->S.prim[VV1] = v*cos(angle);
        c->S.prim[VV2] = v*sin(angle);
        c->S.prim[TR1] = 1;
      }
      else{
        c->S.prim[RHO] = rho_ref;
        c->S.prim[PPP] = p_ref;
        c->S.prim[VV1] = v_ref*cos(angle);
        c->S.prim[VV2] = v_ref*sin(angle);
        c->S.prim[TR1] = 0;
      }
    }

  }

  return 0;

}


void Grid::userKinematics(int it, double t){

  UNUSED(it);
  UNUSED(t);

}

void Cell::userSourceTerms(double dt){

  UNUSED(dt);

}

void Grid::userBoundaries(int it, double t){

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  UNUSED(j);
  UNUSED(i);
  return(skip_);
  
}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value
  
  UNUSED(*res);

}


void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  UNUSED(*p);

  return;

}



void Simu::dataDump(){

  // if (it%1 == 0){ grid.printCols(it, t); }
  if (it%100 == 0){ grid.printCols(it, t); }

}

void Simu::runInfo(){

  // if ((worldrank == 0) and (it%1 == 0)){ printf("it: %ld time: %le\n", it, t);}
  if ((worldrank == 0) and (it%100 == 0)){ printf("it: %ld time: %le\n", it, t);}

}

void Simu::evalEnd(){

  // if (it > 300){ stop = true; }
  if (t > 0.7){ grid.printCols(it, t); stop = true; } // 3.33e8 BOXFIT simu

}



