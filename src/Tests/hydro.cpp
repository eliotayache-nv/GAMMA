/*
* @Author: eliotayache
* @Date:   2020-05-12 10:49:50
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-07-24 11:35:38
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <gsl/gsl_roots.h>   // root finding algorithms
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <mpi.h>

#define X_  1
#define Y_  2

#define HLL_   1
#define HLLC_  2

#define KH_    1
#define TDRP_  2

#define OUTFLOW_  1
#define PERIODIC_ 2

#define ENABLED_   0
#define DISABLED_  1

// ---------------------------------------------------------------------------------------

#define GAMMA_ (5./3.)
#define Nx_    10    // power of 2
#define Ny_    10    // power of 2
#define Ng_    1      // no. of ghost tracks on each side
#define itmax_ 1
#define cfl_   0.2

#define SOLVER_   HLLC_ 
#define CASE_     TDRP_
#define BOUNDARY_ OUTFLOW_

#define OPEN_MPI_ ENABLED_


// ---------------------------------------------------------------------------------------
// Global variables
int numprocs, myid, nodesize, noderank;     // MPI variables
MPI_Comm nodecom;
int nodeNy, nodeNy_gt;

// ---------------------------------------------------------------------------------------
using namespace std;

class Interface;

template <class T> T** array_2d(const int ni, const int nj)
{
  T** pp;

  try {
  pp = new T*[ni];
  }
  catch (bad_alloc&) {
  cout<<"arr_2d: memory allocation error, pp = new T*["<<ni<<"]"<<endl;
  exit(1);
  }

  try {
  pp[0] = new T[ni*nj];
  }
  catch (bad_alloc&) {
  cout<<"arr_2d: memory allocation error, pp[0] = new T["<<ni*nj<<"]"<<endl;
  exit(1);
  }

  for (int i=1; i<ni; ++i) { pp[i] = pp[0] + i * nj; }

  return pp;
}

template <class T> T** to_array_2d(T* A, const int ni, const int nj)
{
  T** pp;

  try {
  pp = new T*[ni];
  }
  catch (bad_alloc&) {
  cout<<"arr_2d: memory allocation error, pp = new T*["<<ni<<"]"<<endl;
  exit(1);
  }

  pp[0] = A;

  for (int i=1; i<ni; ++i) { pp[i] = pp[0] + i * nj; }

  return pp;
}

template <class T> void delete_array_2d(T** pp)
{
  delete [] pp[0];
  delete [] pp;
}

double v2lfac(double v)
{
  return(1. / sqrt(1. - v * v));
}

class Cell
{
public:
  Cell(){}
  ~Cell(){}

  double x, y, rho, vx, vy, p, h, lfac, v, D, mx, my, m, E;   
  double FxD, Fxmx, Fxmy, FxE;
  double FyD, Fymx, Fymy, FyE;
  double dV;

  void update(Interface IL, Interface IR, double dt);

  void prim2cons()
  {
    v = sqrt(vx*vx + vy*vy);
    lfac = v2lfac(v);
    D = rho * lfac;
    h = 1. + p / rho * GAMMA_ / (GAMMA_-1.);
    mx = D * h * lfac * vx;
    my = D * h * lfac * vy;
    E = D * h * lfac - p;
    m = sqrt(mx*mx + my*my);
  }

  void state2flux()
  {
    FxD  = D * vx;
    Fxmx = mx * vx + p;
    Fxmy = my * vx;
    FxE  = mx;

    FyD  = D * vy;
    Fymx = mx * vy;
    Fymy = my * vy + p;
    FyE  = my;
    return;
  }


  struct f_params
  {
    double D,m,E,gamma;
  };


  // This is the function we need to set to zero:
  static double f(double p, void *params)
  {

    // variables
    double lfac;        // Lorentz factor

    // parsing parameters
    struct f_params *par = (struct f_params *) params;
    double  D = par->D;
    double  m = par->m;
    double  E = par->E;
    double  gamma = par->gamma;

    // operations
    lfac = 1. / sqrt(1. - (m * m) / ((E + p) * (E + p)));

    return (E + p - D * lfac - (gamma * p * lfac * lfac) / (gamma - 1.));
      // Mignone (2006) eq. 5
  }

  void cons2prim(double pin = 0)
  {
    double  lfac;                   // Lorentz factor (computed from conserved variables)
    int     status;
    int     iter = 0, max_iter = 1000000;
    double  r;
    double  p_lo, p_hi;
    struct  f_params                params;
    const   gsl_root_fsolver_type   *T;
    gsl_root_fsolver                *s;
    gsl_function                    F;

    F.function = &f;
    F.params = &params;

    // setting initial parameters

    m = sqrt(mx*mx + my*my);

    if (pin != 0) p = pin;
    params.D = D;
    params.m = m;
    params.E = E;
    params.gamma = GAMMA_;

    // Looking for pressure only if current one doesn't work;
    double f_init = f(p, &params);

    if (fabs(f_init) > 1.e-14)
    {
      T = gsl_root_fsolver_brent;
      s = gsl_root_fsolver_alloc (T);

      // setting boundaries
      // f(p_lo) has to be positive. No solution otherwise
      // hence, f(p_hi) has to be negative

      p_lo = fmax(0., (1. + 1e-13) * fabs(m) - E);
      // printf("%le %le %le %le\n", p_lo, m, fabs(m) - E, f(p_lo, &params));
      if (f_init < 0)
      {
        p_hi = 1. * p;
      }
      else
      {
        int i = 0;

        p_hi = p;
        while (f(p_hi, &params) > 0)
        {
          i++;
          p_hi *= 10;
        }
      }
      gsl_root_fsolver_set (s, &F, p_lo, p_hi);

      // Iterating for root finding:            
      do
        {
          iter++;
          status = gsl_root_fsolver_iterate (s);
          r = gsl_root_fsolver_root (s);
          p_lo = gsl_root_fsolver_x_lower (s);
          p_hi = gsl_root_fsolver_x_upper (s);
          status = gsl_root_test_interval (p_lo, p_hi,
                           0, 1.e-13);
        }                
        while (status == GSL_CONTINUE && iter < max_iter);

      gsl_root_fsolver_free (s);
      p = r;
    }

    lfac = 1. / sqrt(1. - (m * m) / ((E + p) * (E + p)));
    rho = D / lfac;
    if (m == 0.)    
    {
      v  = 0.;
      vx = 0.;
      vy = 0.;
    }
    else
    {
      v   = sqrt(1. - 1. / (lfac * lfac));
      // Mignone (2006) eq. 3
      vx = v * (mx/m);
      vy = v * (my/m);            
    }

    prim2cons(); // for consistency

    return;
  }
};

void generate_mpi_cell( MPI_Datatype * cell_mpi )
{
  Cell cell;
  int count = 1; // no. of types: only doubles
  int numVariables = 23;  // number of variables to communicate (state + volume)
  MPI_Datatype types[count];
  int blocklenghts[count];
  MPI_Aint offsets[count];

  types[0] = MPI_DOUBLE;
  blocklenghts[0] = numVariables;

  MPI_Aint address1, address2;
  MPI_Get_address(&cell,&address1);   // address of cell
  MPI_Get_address(&cell.x,&address2); // address of first attribute of cell

  offsets[0] = address2 - address1;

  MPI_Type_create_struct(count, blocklenghts, offsets, types, cell_mpi);
  MPI_Type_commit(cell_mpi);
}

class Interface
{
public:
  Interface(){}
  ~Interface(){}
  
  double FD, Fmx, Fmy, FE, lL, lR, lStar;
  double dA;

  void computeLambda(Cell cL, Cell cR, int dim = X_)
  {
    double  sigmaSL, sigmaSR;       // speed of sound param
    double  cSL, cSR;               // speed of sound
    double  lambda1, lambda2;       // 2 intermediates
    double  vL, vR, mL, mR, pL, pR, EL, ER;

    if (dim == X_)
    {
      vL = cL.vx;
      vR = cR.vx;
      mL = cL.mx;
      mR = cR.mx;
    }
    else if (dim == Y_)
    {
      vL = cL.vy;
      vR = cR.vy;
      mL = cL.my;
      mR = cR.my;
    }

    pL = cL.p;
    pR = cR.p;
    EL = cL.E;
    ER = cR.E;

    // LEFT AND RIGHT WAVE SPEEDS ESTIMATsES
    // Computing speed of sound parameter sigmaS:
    cSL = sqrt(GAMMA_ * cL.p / (cL.rho + cL.p * GAMMA_ / (GAMMA_ - 1.)));
    cSR = sqrt(GAMMA_ * cR.p / (cR.rho + cR.p * GAMMA_ / (GAMMA_ - 1.)));

      // Mignone (2006) eq. 4
    sigmaSL = cSL * cSL / (cL.lfac * cL.lfac * (1. - cSL * cSL));
    sigmaSR = cSR * cSR / (cR.lfac * cR.lfac * (1. - cSR * cSR));
      // Mignone (2006) eq. 22-23s

    // Retreiving L and R speeds:
    lambda1 = (vR - sqrt(sigmaSR * (1. - vR * vR + sigmaSR))) / (1. + sigmaSR);
    lambda2 = (vL - sqrt(sigmaSL * (1. - vL * vL + sigmaSL))) / (1. + sigmaSL);

    lL = 1. * fmin(lambda1,lambda2);

    lambda1 = (vR + sqrt(sigmaSR * (1. - vR * vR + sigmaSR))) / (1. + sigmaSR);
    lambda2 = (vL + sqrt(sigmaSL * (1. - vL * vL + sigmaSL))) / (1. + sigmaSL);

    lR = 1. * fmax(lambda1,lambda2);

    // for HLLC
    double  AL,AR;                      // temporary variables
    double  BL,BR;                      // temporary variables
    double  FhllE, Ehll, Fhllm, mhll;   // Coefficients of the equation for lambdaStar
      // Mignone (2006) eq. 18
    double delta;                       // equation discriminant

    // Setting up temporary variables:
    AL = lL * EL - mL;
    AR = lR * ER - mR;
    BL = mL * (lL - vL) - pL;
    BR = mR * (lR - vR) - pR;
      // Mignone (2006) eq. 17

    // Coefficients for lambdaStar equation:
    FhllE = (lL * AR - lR * AL) / (lR - lL);
    Ehll  = (AR - AL) / (lR - lL);
    Fhllm = (lL * BR - lR * BL) / (lR - lL);
    mhll  = (BR - BL) / (lR - lL);

    // Solving for lambdaStar:
    // Mignone (2006) eq. 18
    /* ----------------------------------------------------- */
    /* lambdaStar is the (-) solution of the quadratic eq.   */
    /* This is for physical reasons.                         */
    /*                                                       */
    /* The result is tested to be inside the lL,lR interval  */
    /* ----------------------------------------------------- */
    if (FhllE == 0)
    {
      lStar = mhll / (Ehll + Fhllm);
      return;
    }

    delta = pow(Ehll + Fhllm,2.) - 4. * FhllE * mhll;
    lStar = ((Ehll + Fhllm) - sqrt(delta)) / (2. * FhllE);
  }

  void computeF(Cell cL, Cell cR, int dim = X_)
  {
    Cell   cstar;              // hll star state

    if (dim == X_)
    {
      #if SOLVER_ == HLLC_
      //variables:
      double  AL,AR;                      // temporary variables
      double  BL,BR;                      // temporary variables
      Cell    cLstar, cRstar;     // hllc star states        

      // operations:
      // Computing star region states

      // Setting up temporary variables:
      AL = lL * cL.E - cL.mx;
      AR = lR * cR.E - cR.mx;
      BL = cL.mx * (lL - cL.vx) - cL.p;
      BR = cR.mx * (lR - cR.vx) - cR.p;
        // Mignone (2006) eq. 17

      cLstar.p = (AL * lStar - BL) / (1. - lL * lStar);
      cRstar.p = (AR * lStar - BR) / (1. - lR * lStar);
        // Mignone (2006) eq. 17

      // if (lL != 0 || lR != 0 || lStar != 0)
      // U*left
      cLstar.D = cL.D * (lL - cL.vx) / (lL - lStar);
      cLstar.mx = (cL.mx * (lL - cL.vx) + cLstar.p - cL.p) / (lL - lStar); 
      cLstar.my = cL.my * (lL - cL.vx) / (lL - lStar);
      cLstar.E = (cL.E * (lL - cL.vx) + cLstar.p * lStar - cL.p * cL.vx) 
        / (lL - lStar);
          // Mignone (2006) eq. 16

      // U*right
      cRstar.D = cR.D * (lR - cR.vx) / (lR - lStar);
      cRstar.mx = (cR.mx * (lR - cR.vx) + cRstar.p - cR.p) / (lR - lStar);
      cRstar.my = cR.my * (lR - cR.vx) / (lR - lStar);
      cRstar.E = (cR.E * (lR - cR.vx) + cRstar.p * lStar - cR.p * cR.vx) 
        / (lR - lStar);
          // Mignone (2006) eq. 16

      // Computing star states primitive velocity:
      cLstar.vx = lStar;
      cRstar.vx = lStar;

      // Computing fluxes
      cLstar.state2flux();
      cRstar.state2flux();

      // Choosing right flux:
      // when solving in the interface rest frame
      if (0 < lL)
      {
        FD = cL.FxD;
        Fmx = cL.Fxmx;
        Fmy = cL.Fxmy;
        FE = cL.FxE;
      }                     
      if (lL <= 0 && 0 < lStar)
      {   
        FD = cLstar.FxD;
        Fmx = cLstar.Fxmx;
        Fmy = cLstar.Fxmy;
        FE = cLstar.FxE;
      }
      if (lStar <= 0 && 0 < lR)
      {
        FD = cRstar.FxD;
        Fmx = cRstar.Fxmx;
        Fmy = cRstar.Fxmy;
        FE = cRstar.FxE;
      }
      if (lR <= 0 )
      {
        FD = cR.FxD;
        Fmx = cR.Fxmx;
        Fmy = cR.Fxmy;
        FE = cR.FxE;
      }
      #endif

      #if SOLVER_ == HLL_
      // Computing HLL region state
      cstar.D  = (lR * cR.D - lL * cL.D + cL.FxD - cR.FxD) / (lR - lL);
      cstar.mx = (lR * cR.mx - lL * cL.mx + cL.Fxmx - cR.Fxmx) / (lR - lL);
      cstar.my = (lR * cR.my - lL * cL.my + cL.Fxmy - cR.Fxmy) / (lR - lL);
      cstar.E  = (lR * cR.E - lL * cL.E + cL.FxE - cR.FxE) / (lR - lL);
        // Mignone (2006) eq. 9

      // Computing HLL region flux
      cstar.FxD  = (lR * cL.FxD - lL * cR.FxD + lR * lL * (cR.D - cL.D)) / (lR - lL);
      cstar.Fxmx = (lR * cL.Fxmx - lL * cR.Fxmx + lR * lL * (cR.mx - cL.mx)) / (lR - lL);
      cstar.Fxmy = (lR * cL.Fxmy - lL * cR.Fxmy + lR * lL * (cR.my - cL.my)) / (lR - lL);
      cstar.FxE  = (lR * cL.FxE - lL * cR.FxE + lR * lL * (cR.E - cL.E)) / (lR - lL);
        // Mignone (2006) eq. 11

      // Choosing right flux:
      // when solving in the interface rest frame
      if (0 < lL)
      {
        FD = cL.FxD;
        Fmx = cL.Fxmx;
        Fmy = cL.Fxmy;
        FE = cL.FxE;
      }                     
      if (lL <= 0 && 0 < lR)
      {   
        FD = cstar.FxD;
        Fmx = cstar.Fxmx;
        Fmy = cstar.Fxmy;
        FE = cstar.FxE;
      }
      if (lR <= 0 )
      {
        FD = cR.FxD;
        Fmx = cR.Fxmx;
        Fmy = cR.Fxmy;
        FE = cR.FxE;
      }
      #endif
    }
    else if (dim == Y_)
    {

      #if SOLVER_ == HLLC_
      //variables:
      double  AL,AR;                      // temporary variables
      double  BL,BR;                      // temporary variables
      Cell    cLstar, cRstar;     // hllc star states        

      // operations:
      // Computing star region states

      // Setting up temporary variables:
      AL = lL * cL.E - cL.my;
      AR = lR * cR.E - cR.my;
      BL = cL.my * (lL - cL.vy) - cL.p;
      BR = cR.my * (lR - cR.vy) - cR.p;
        // Mignone (2006) eq. 17

      cLstar.p = (AL * lStar - BL) / (1. - lL * lStar);
      cRstar.p = (AR * lStar - BR) / (1. - lR * lStar);
        // Mignone (2006) eq. 17

      // if (lL != 0 || lR != 0 || lStar != 0)
      // U*left
      cLstar.D = cL.D * (lL - cL.vy) / (lL - lStar);
      cLstar.my = (cL.my * (lL - cL.vy) + cLstar.p - cL.p) / (lL - lStar); 
      cLstar.mx = cL.mx * (lL - cL.vy) / (lL - lStar); 
      cLstar.E = (cL.E * (lL - cL.vy) + cLstar.p * lStar - cL.p * cL.vy) 
        / (lL - lStar);
          // Mignone (2006) eq. 16

      // U*right
      cRstar.D = cR.D * (lR - cR.vy) / (lR - lStar);
      cRstar.my = (cR.my * (lR - cR.vy) + cRstar.p - cR.p) / (lR - lStar);
      cRstar.mx = cR.mx * (lR - cR.vy) / (lR - lStar);
      cRstar.E = (cR.E * (lR - cR.vy) + cRstar.p * lStar - cR.p * cR.vy) 
        / (lR - lStar);
          // Mignone (2006) eq. 16

      // Computing star states primitive velocity:
      cLstar.vy = lStar;
      cRstar.vy = lStar;

      // Computing fluxes
      cLstar.state2flux();
      cRstar.state2flux();

      // Choosing right flux:
      // when solving in the interface rest frame
      if (0 < lL)
      {
        FD = cL.FyD;
        Fmx = cL.Fymx;
        Fmy = cL.Fymy;
        FE = cL.FyE;
      }                     
      if (lL <= 0 && 0 < lStar)
      {   
        FD = cLstar.FyD;
        Fmx = cLstar.Fymx;
        Fmy = cLstar.Fymy;
        FE = cLstar.FyE;
      }
      if (lStar <= 0 && 0 < lR)
      {
        FD = cRstar.FyD;
        Fmx = cRstar.Fymx;
        Fmy = cRstar.Fymy;
        FE = cRstar.FyE;
      }
      if (lR <= 0 )
      {
        FD = cR.FyD;
        Fmx = cR.Fymx;
        Fmy = cR.Fymy;
        FE = cR.FyE;
      }
      #endif

      #if SOLVER_ == HLL_
      // Computing HLL region state
      cstar.D  = (lR * cR.D - lL * cL.D + cL.FyD - cR.FyD) / (lR - lL);
      cstar.mx = (lR * cR.mx - lL * cL.mx + cL.Fymx - cR.Fymx) / (lR - lL);
      cstar.my = (lR * cR.my - lL * cL.my + cL.Fymy - cR.Fymy) / (lR - lL);
      cstar.E  = (lR * cR.E - lL * cL.E + cL.FyE - cR.FyE) / (lR - lL);
        // Mignone (2006) eq. 9

      // Computing HLL region flux
      cstar.FyD  = (lR * cL.FyD - lL * cR.FyD + lR * lL * (cR.D - cL.D)) / (lR - lL);
      cstar.Fymx = (lR * cL.Fymx - lL * cR.Fymx + lR * lL * (cR.mx - cL.mx)) / (lR - lL);
      cstar.Fymy = (lR * cL.Fymy - lL * cR.Fymy + lR * lL * (cR.my - cL.my)) / (lR - lL);
      cstar.FyE  = (lR * cL.FyE - lL * cR.FyE + lR * lL * (cR.E - cL.E)) / (lR - lL);
        // Mignone (2006) eq. 11            

      // Choosing right flux:
      // when solving in the interface rest frame
      if (0 < lL)
      {
        FD = cL.FyD;
        Fmx = cL.Fymx;
        Fmy = cL.Fymy;
        FE = cL.FyE;
      }                     
      if (lL <= 0 && 0 < lR)
      {   
        FD = cstar.FyD;
        Fmx = cstar.Fymx;
        Fmy = cstar.Fymy;
        FE = cstar.FyE;
      }
      if (lR <= 0 )
      {
        FD = cR.FyD;
        Fmx = cR.Fymx;
        Fmy = cR.Fymy;
        FE = cR.FyE;
      }
      #endif
    }
  }
};

void Cell::update(Interface IL, Interface IR, double dt)
{
  D += (IL.FD*IL.dA - IR.FD*IR.dA) * dt / dV;
  mx += (IL.Fmx*IL.dA - IR.Fmx*IR.dA) * dt / dV;
  my += (IL.Fmy*IL.dA - IR.Fmy*IR.dA) * dt / dV;
  E += (IL.FE*IL.dA - IR.FE*IR.dA) * dt / dV;
}

void print2Darray(Cell **C, const int nx = Nx_, const int ny = Ny_)
{
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      printf("%le ", C[j][i].rho);
    }
    printf("\n");
  }
}

void save2Darray(Cell **C, const int it, const int nx = Nx_, const int ny = Ny_)
{
  ofstream frho, fvx, fvy, fp;
  ofstream fD, fmx, fmy, fE;
  string filename_rho, filename_vx, filename_vy, filename_p;
  string filename_D, filename_mx, filename_my, filename_E;

  filename_rho = "data/rho" + to_string(it) + ".out";
  filename_vx = "data/vx" + to_string(it) + ".out";
  filename_vy = "data/vy" + to_string(it) + ".out";
  filename_p = "data/p" + to_string(it) + ".out";
  filename_D = "data/D" + to_string(it) + ".out";
  filename_mx = "data/mx" + to_string(it) + ".out";
  filename_my = "data/my" + to_string(it) + ".out";
  filename_E = "data/E" + to_string(it) + ".out";

  frho.open(filename_rho, ios::trunc);
  fvx.open(filename_vx, ios::trunc);
  fvy.open(filename_vy, ios::trunc);
  fp.open(filename_p, ios::trunc);
  fD.open(filename_D, ios::trunc);
  fmx.open(filename_mx, ios::trunc);
  fmy.open(filename_my, ios::trunc);
  fE.open(filename_E, ios::trunc);

  frho << scientific;
  fvx << scientific;
  fvy << scientific;
  fp << scientific;
  fD << scientific;
  fmx << scientific;
  fmy << scientific;
  fE << scientific;

  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      frho << C[j][i].rho << " ";
      fvx << C[j][i].vx << " ";
      fvy << C[j][i].vy << " ";
      fp << C[j][i].p << " ";
      fD << C[j][i].D << " ";
      fmx << C[j][i].mx << " ";
      fmy << C[j][i].my << " ";
      fE << C[j][i].E << " ";
    }
    frho << endl;
    fvx << endl;
    fvy << endl;
    fp << endl;
    fD << endl;
    fmx << endl;
    fmy << endl;
    fE << endl;
  }
  frho.close();
  fvx.close();
  fvy.close();
  fp.close();
  fD.close();
  fmx.close();
  fmy.close();
  fE.close();
}

double timescale(Interface **I, int dir=X_)
{
  double dt = 1.e15;
  int maxj;
  if (dir==X_) maxj=nodeNy; else if (dir==Y_) maxj=nodeNy+1;
  for (int j = 0; j < maxj; ++j)
    for (int i = 0; i < Nx_; ++i)
    {
      dt = fmin(dt, fmin(fabs(1./(Nx_*I[j][i].lL)), fabs(1./(Nx_*I[j][i].lR))));
    }
  return cfl_ * dt;
}


void fillGhostTracks(Cell **Cgrid, Cell **C)
{
  MPI_Datatype cell_mpi = {0}; 
  generate_mpi_cell( &cell_mpi );

  int below_id = (myid + numprocs - nodesize) % numprocs;  // ring communication
  //                     ^ added to shift everyone to postitive values
  int above_id = (myid + nodesize) % numprocs;  // overidden by boundary conditions
  int node_num = myid / nodesize;
  int node_count = numprocs / nodesize;

  MPI_Barrier(MPI_COMM_WORLD);
  if (noderank == 0) 
  {
    // lower ghost track comes from lower node
    MPI_Sendrecv(&C[nodeNy-1][0], Nx_, cell_mpi, above_id, 0, 
                 &Cgrid[0][0]   , Nx_, cell_mpi, below_id, 0, 
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // higher ghost track comes from higher node
    MPI_Sendrecv(&C[0][0]              , Nx_, cell_mpi, below_id, 1, 
                 &Cgrid[nodeNy_gt-1][0], Nx_, cell_mpi, above_id, 1, 
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


int main(int argc, char *argv[])
{
  double t, dt, local_dt;
  Cell *Cmem, **C, **Cgrid;
  Interface *Ixmem, **Ix;     // interfaces in the x direction
  Interface *Iymem, **Iy;     // interfaces in the y direction

  // Initialising MPI if needed
  #if OPEN_MPI_ == ENABLED_

    MPI_Aint Csize, Ixsize, Iysize;
    MPI_Aint Clocalsize, Ixlocalsize, Iylocalsize; 
      // only noderank 0 gets to allocate the full array
    int Csize_disp, Ixsize_disp, Iysize_disp;

    // dynamic RMA windows for shared memory allocation
    MPI_Win Cwin, Ixwin, Iywin;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up

    // get node-specific information
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, myid,
      MPI_INFO_NULL, &nodecom);

    MPI_Comm_size(nodecom, &nodesize);
    MPI_Comm_rank(nodecom, &noderank);
  #elif OPEN_MPI_ == DISABLED_
    numprocs = 1;
    myid     = 0;
    noderank = 0;
  #endif

  if (myid == 0)
  {
    printf("\nMPI STATUS:\n");
    printf("numprocs: %d\n", numprocs);
    printf("nodesize: %d\n\n", nodesize);
  }


  #if OPEN_MPI_ == ENABLED_

    nodeNy = Ny_ / (numprocs / nodesize); // we're splitting in y-direction
    nodeNy_gt = nodeNy + 2*Ng_; // with ghost tracks

    if (noderank == 0)
    {
      Csize = Nx_ * nodeNy_gt;
      Ixsize = (Nx_+1) * nodeNy;
      Iysize = Nx_ * (nodeNy+1);
      Clocalsize = Csize;
      Ixlocalsize = Ixsize;
      Iylocalsize = Iysize;
    }
    else
    {
      Clocalsize = 0;
      Ixlocalsize = 0;
      Iylocalsize = 0;      
    }

    MPI_Win_allocate_shared(Clocalsize*sizeof(Cell), sizeof(Cell),
      MPI_INFO_NULL, nodecom, &Cmem, &Cwin);
    MPI_Win_allocate_shared(Ixlocalsize*sizeof(Interface), sizeof(Interface),
      MPI_INFO_NULL, nodecom, &Ixmem, &Ixwin);
    MPI_Win_allocate_shared(Iylocalsize*sizeof(Interface), sizeof(Interface),
      MPI_INFO_NULL, nodecom, &Iymem, &Iywin);

    // tell cores not of top rank in their respective node where to find data
    if (noderank != 0)
    {    
      MPI_Win_shared_query(Cwin, 0, &Csize, &Csize_disp, &Cmem);
      MPI_Win_shared_query(Ixwin, 0, &Ixsize, &Ixsize_disp, &Ixmem);
      MPI_Win_shared_query(Iywin, 0, &Iysize, &Iysize_disp, &Iymem);
    }

    MPI_Barrier(nodecom); // wait until all cores have caught up

    Cgrid = to_array_2d<Cell>(Cmem, nodeNy_gt, Nx_);   // includes ghost tracks
    C  = &Cgrid[Ng_];  // offset by no. of ghost tracks (only includes physical grid)
    Ix = to_array_2d<Interface>(Ixmem, nodeNy, Nx_+1);
    Iy = to_array_2d<Interface>(Iymem, nodeNy+1, Nx_);
      // The "track" dimension has to be the deepest
      // i is iter for x
      // j ----------- y
      // Ix are interfaces between cells of same x
      // Iy ------------------------------------ y
  #else
    // Broken by ghost tracks update
    nodeNy = Ny_;
    C  = array_2d<Cell>(Ny_, Nx_);   // wrapper to access big mem location in 2D
    Ix = array_2d<Interface>(Ny_, Nx_+1);
    Iy = array_2d<Interface>(Ny_+1, Nx_);
  #endif

  // Initial conditions
  if (noderank == 0)
  {
    int basej = myid / nodesize * nodeNy; // again, only works with grids as powers of 2
      // or at least nodeNy is the same in all nodes
    for (int j = 0; j < nodeNy; ++j)
    {
      for (int i = 0; i < Nx_; ++i)
      {
        Cell *pc = &C[j][i];

        pc->x = (double) 2. * (i + 0.5) / Nx_ - 1.;
        pc->y = (double) 2. * (basej + j + 0.5) / Ny_ - 1.; // basej for each node
        pc->dV = 4. / (Nx_*Ny_);

        #if CASE_ == KH_

          pc->vy = 0.;

          if (pc->y > 0.) 
          {
            pc->vx = 0.1;
            pc->rho = 2.;
            pc->p = 2.5;
          }
          else 
          {
            pc->vx = -0.1;
            pc->rho = 1.;
            pc->p = 2.5;
          }

          // pc->vx += 0.01 * (double) gsl_rng_get (r) / (double) gsl_rng_max(r);
          // pc->vy += 0.01 * (double) gsl_rng_get (r) / (double) gsl_rng_max(r);

        #elif CASE_ == TDRP_

          if (pc->x > 0 and pc->y > 0)
          {
            pc->rho = 0.1;
            pc->vx  = 0.;
            pc->vy  = 0.;
            pc->p   = 0.01;
          }
          if (pc->x < 0 and pc->y > 0)
          {
            pc->rho = 0.1;
            pc->vx  = 0.99;
            pc->vy  = 0.;
            pc->p   = 1.;
          }
          if (pc->x < 0 and pc->y < 0)
          {
            pc->rho = 0.5;
            pc->vx  = 0.;
            pc->vy  = 0.;
            pc->p   = 1.;
          }
          if (pc->x > 0 and pc->y < 0)
          {
            pc->rho = 0.1;
            pc->vx  = 0.;
            pc->vy  = 0.99;
            pc->p   = 1.;
          }

        #endif

        pc->prim2cons();
        pc->state2flux();
      }
    }
  }

  // Splitting computation zones based on noderank
  int jmin = noderank * nodeNy / nodesize;
  int jmax = (noderank+1) * nodeNy / nodesize -1;
  int jmax_I = jmax;
  if (noderank == nodesize-1) jmax_I++; // extra boundary interface 

  // Evolution
  t = 0;
  for (int it = 0; it < itmax_; ++it)
  {
    MPI_Barrier(MPI_COMM_WORLD);

    // Writing out
    if (it%20 == 0)
    { 
      int size  = Ny_ / numprocs * Nx_;
      MPI_Datatype cell_mpi = {0}; 
      generate_mpi_cell( &cell_mpi );

      if (myid == 0)
      {
        Cell **Cdump  = array_2d<Cell>(Ny_, Nx_);   // wrapper to access big mem location in 2D

        // copying info from proc 0
        std::memcpy(&Cdump[0][0],&C[0][0],size * sizeof(Cell));

        // receiving from all non-0 processes
        for (int j = 1; j < numprocs; ++j)
        {
          int index = j * Ny_ / numprocs;
          MPI_Recv(&Cdump[index][0], size, cell_mpi, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        cout << t << endl;
        save2Darray(Cdump, it);
        delete_array_2d<Cell>(Cdump);        
      } 
      else
      {
        int index = noderank * nodeNy / nodesize;
        MPI_Send(&C[index][0], size, cell_mpi, 0, 0, MPI_COMM_WORLD);
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }

    fillGhostTracks(Cgrid, C);

    for (int j = jmin; j <= jmax_I; ++j)
    {
      for (int i = 0; i <= Nx_; ++i)
      {
        // Linking neighbours
        int iL, iR, jL, jR;

        #if BOUNDARY_ == PERIODIC_

          if (i == Nx_ or j == Ny_) continue;

          iR = i;
          if (i!=0)   iL = i - 1;
          else        iL = Nx_-1;
          jR = j;
          if (j!=0)   jL = j - 1;
          else        jL = 0;

        #elif BOUNDARY_ == OUTFLOW_

          if (i!=Nx_) iR = i;
          else        iR = i-1;
          if (i!=0)   iL = i - 1;
          else        iL = 0;
          // if (j!=Ny_) jR = j;
          // else        jR = j-1;
          // if (j!=0)   jL = j - 1;
          // else        jL = 0;

          // copying tracks for outflow boundary condition
          if (myid==0) 
            std::memcpy(&Cgrid[0][0],&C[0][0],Nx_ * sizeof(Cell));
          if (myid==numprocs-1)
            std::memcpy(&Cgrid[nodeNy_gt-1][0],&C[nodeNy-1][0],Nx_ * sizeof(Cell));

        #endif

        if (j!=nodeNy)
        {
          Ix[j][i].dA = 2./Ny_;
          Ix[j][i].computeLambda(C[j][iL], C[j][iR], X_);
          Ix[j][i].computeF(C[j][iL], C[j][iR], X_);
        }

        if (i!=Nx_)
        {
          Iy[j][i].dA = 2./Nx_;
          Iy[j][i].computeLambda(C[j-1][i], C[j][i], Y_);
          Iy[j][i].computeF(C[j-1][i], C[j][i], Y_);
        }
      }
    }


    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up

    // computing dt
    local_dt = fmin(timescale(Ix, X_), timescale(Iy, Y_));
    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
    // taking minimum dt accross all processes
    MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    t += dt;

    for (int j = jmin; j <= jmax; ++j)
    {

      for (int i = 0; i < Nx_; ++i)
      {
        // Linking neighbours
        int iL, iR, jL, jR;

        #if BOUNDARY_ == PERIODIC_

          if (i == Nx_ or j == Ny_) continue;

          iL = i;
          if (i!=Nx_-1) iR = i+1;
          else          iR = 0;
          jL = j;
          if (j!=Ny_-1) jR = j+1;
          else          jR = Ny_-1;

        #elif BOUNDARY_ == OUTFLOW_

          iL = i;
          iR = i+1;
          jL = j;
          jR = j+1;

        #endif

        C[j][i].update(Ix[j][iL], Ix[j][iR], dt);
        C[j][i].update(Iy[jL][i], Iy[jR][i], dt);

        C[j][i].cons2prim();
        C[j][i].state2flux();
      }
    }
  }

  #if OPEN_MPI_ == ENABLED_

    MPI_Barrier(MPI_COMM_WORLD); // wait until all cores have caught up
    MPI_Finalize();

  #endif // OPEN_MPI_ == ENABLED_

  return 0;
}