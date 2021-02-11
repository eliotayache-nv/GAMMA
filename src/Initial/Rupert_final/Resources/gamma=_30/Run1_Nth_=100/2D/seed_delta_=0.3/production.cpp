#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

// initial Fireball (FB) setup

static double lfac    = 30.;                 // --:      energy per unit mass ~ coasting Lorentz factor
static double Eiso    = 1.e52;                // erg:     isotropic equivalent energy
static double R0      = 100.;                  // ls:      initial shell width
static double n0      = 1.;                   // cm-3:    CBM proton number density

// simulation settings
static double th_simu = PI/32.;        // rad:     simulation angle
static double t_ini   = 0.;                   // s:     start time of dynamic simulation
static double rmin    = 0.01;                  // light-seconds, box inner boundary at t_ini
static double rmax    = 150.;                 // light-seconds, box outer boundary at t_ini

static double R03     = R0*R0*R0;
static double c2      = c_*c_;
static double c3      = c_*c_*c_;

// set external medium density-pressure relation
static double eta     = 1.e-5;                // p0 = rho0*c2, since cold external medium

// calculate CBM parameters
static double rho0    = n0*mp_;               // g.cm-3:  CBM lab frame mass density
static double p0      = rho0*eta*c2;

// calculate FB parameters
static double Miso    = Eiso/(lfac*c2);       // g:         FB initial isotropic equivalent mass
static double Viso    = (4./3.)*PI*R03*c3;    // cm3:       FB initial isotropic equivalent volume
static double E_f     = Eiso/Viso;            // erg.cm-3:  FB total energy density: E = tau + D
static double D_f     = Miso/Viso;            // g.cm-3:    FB mass density (lab frame, but initial velocity is 0, so same in comoving frame)
static double E_D_f   = D_f * c2;             // erg.cm-3:  FB rest-mass energy density
static double tau_f   = E_f - E_D_f;          // erg.cm-3:  FB internal energy density

// normalisation constants:
static double rhoNorm   = rho0;               // density normalised to CBM density
static double lNorm   = c_;                   // distance normalised to c
static double vNorm   = c_;                   // distance normalised to c
static double pNorm   = rhoNorm*vNorm*vNorm;  // pressure normalised to rho_0*c2
static double ENorm   = rhoNorm*lNorm*lNorm;  // energy density normalised based on E/V = E/L3 = (M/L3)*L2/T2 = D*L2/T2, T = time, no normalisation

// regridding parameters
static double fin_AR      = 1.;               //      final target AR
static double R1          = 1.e4;             // ls:  radius to use for target_AR calculation: target_AR = fmax(pow(r_AR_1/r, 0.2), 1.)*fin_AR;
static double R2          = 2.e6;
static double split_AR    = 5.;               //      set upper bound on AR as ratio of target_AR
static double merge_AR    = 0.2;              //      set upper bound on AR as ratio of target_AR
// static double tMove       = 4.e3;             // s:   time at which inner boundary starts to move

// external medium perturbation
static double delta   = 0.3;

void loadParams(s_par *par){

  par->tini      = t_ini;             // initial time
  par->ncell[x_] = 300;               // number of cells in r direction
  par->ncell[y_] = 100;               // number of cells in theta direction
  par->nmax      = 1500;              // max number of cells in MV direction
  par->ngst      = 2;                 // number of ghost cells

}

double random_number(){							                        // random number generator function
	
    // generates random numbers in range (0,1) using minimal standard generator (uniform probability)
    
  static int64_t X_n = 1; 					                        // declare seed as 1; static means it will be retained for future function calls
	int32_t a,c,m;								                        // declare a,c,m as long ints
	a = 16807;									                        // set value of a
	m = 2147483647;								                        // set value of m
  c = 0;										                        // set value of c
	X_n = (a*X_n+c)%m;							                        // calculate random number
	
	double ud;									                        // declare variable as double
	ud = (double)X_n/(m-1);						                        // calculate uniform deviates (ud) in range (0,1)
    
  return ud;									                        // return ud in (0,1)
}

int Grid::initialGeometry(){                              // loop across grid, calculate cell positions and grid spacing, add info to cell struct "c"

// grid defined in length units of light-seconds

  double logrmin  = log(rmin);
  double logrmax  = log(rmax);
  double dlogr    = (logrmax - logrmin);
  
  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      
      double rlog   = (double) dlogr*(i+0.5)/ncell[x_] + logrmin;
      double rlogL  = (double) dlogr*(i  )/ncell[x_] + logrmin;
      double rlogR  = (double) dlogr*(i+1)/ncell[x_] + logrmin;
      double r_ls   = exp(rlog);
      double dr_ls  = exp(rlogR) - exp(rlogL);
      
      double th     = (double) th_simu*(j+0.5)/ncell[y_];
      double dth    = th_simu/ncell[y_];

      c->G.x[x_]    = r_ls;
      c->G.dx[x_]   = dr_ls;
      c->G.x[y_]    = th;
      c->G.dx[y_]   = dth;

      c->computeAllGeom();
    }
  }
  return 0;
}

int Grid::initialValues(){

  if (GAMMA_ != 4./3.){
    printf("WARNING - Set GAMMA_ to 4./3.\n");
  }
  if (VI != 1.){
    printf("WARNING - Set VI to 1.\n");
  }

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){        // loop through cells along theta
    for (int i = 0; i < ncell[MV]; ++i){      // loop through cells along r
      Cell *c = &Cinit[j][i];

      double r = c->G.x[x_];                  // ls:  radial coordinate
      // double th = c->G.x[y_];                 // rad: theta angular coordinate

      if (r <= R0){
        c->S.cons[DEN] = D_f/rhoNorm;
        c->S.cons[TAU] = tau_f/ENorm;
        c->S.cons[SS1] = 0.;
        c->S.cons[SS2] = 0.;
        c->S.cons[TR1] = 1.;
        c->S.cons2prim(r);
        // printf("rho = %e, p = %e\n",c->S.prim[RHO],c->S.prim[PPP]);
        // printf("rho = %e, p = %e\n",(c->S.prim[RHO])*rhoNorm,(c->S.prim[PPP])*pNorm);
        // double gam = (c->S.cons[TAU]+c->S.cons[DEN])/c->S.cons[DEN];
        // printf("Lorentz Factor = %e\n",gam);
      }
      else{
        c->S.prim[RHO] = rho0/rhoNorm;
        c->S.prim[PPP] = p0/pNorm;
        c->S.prim[VV1] = 0.;
        c->S.prim[VV2] = 0.;
        c->S.prim[TR1] = 2.;
        // printf("rho = %e, p = %e\n",c->S.prim[RHO],c->S.prim[PPP]);
      }
    }
  }

  return 0;
}


void Grid::userKinematics(int it, double t){

  // setting lower and higher i boundary interface velocities

  // set by boundary velocity:
  // double vIn;
  // if (t < tMove){
  //   vIn = 0.;
  // }
  // else{
  //   vIn = 0.55;
  // }

  double vIn = 0.;
  double vOut = 1.05;

  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n <= ngst; ++n){
      int    iL = n;
      int    iR = ntrack[j]-2-n;
      Itot[j][iL].v = vIn;
      Itot[j][iR].v = vOut;
    }
  }
}

void Cell::userSourceTerms(double dt){

  UNUSED(dt);

}

void Grid::userBoundaries(int it, double t){

  // outflow BC at inner boundary

  // fixed value BC at outer boundary
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int i = 0; i < ngst; ++i){

      double perturbation = 2.*delta*(random_number() - 0.5);
      double dens = (rho0/rhoNorm)*(1. + perturbation);

      int nt = ntrack[j];
      // Ctot[j][nt-1-i] : outer boundary ghost cells
      // set outer boundary to fixed values of external medium
      Ctot[j][nt-1-i].S.prim[RHO] = dens;         // outer boundary
      Ctot[j][nt-1-i].S.prim[PPP] = p0/pNorm;             // outer boundary
      Ctot[j][nt-1-i].S.prim[VV1] = 0.;                   // outer boundary
      Ctot[j][nt-1-i].S.prim[VV2] = 0.;                   // outer boundary
      Ctot[j][nt-1-i].S.prim[TR1] = 2.;                   // outer boundary
    }
  }
  
  // reflective BCs on inner side (along jet axis)

  if (worldrank==0){
    for (int j = 0; j < ngst; ++j){

      int target_j = 2*ngst-1-j;
      ntrack[j] = ntrack[target_j];
      nact[j] = nact[target_j];
      iRbnd[j] = iRbnd[target_j];
      iLbnd[j] = iLbnd[target_j];
      std::copy_n(&Ctot[target_j][0], ntrack[target_j],   &Ctot[j][0]);
      std::copy_n(&Itot[target_j][0], ntrack[target_j]-1, &Itot[j][0]);

      // updating ghost positions, ids and indexes
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] -= (jLbnd-j+1)*C[0][0].G.dx[F1]; 
          Itot[j][i].computedA();
        }
      }

      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  // reflective BCs on outer side (outside of jet cone)

  if (worldrank==worldsize-1){
    for (int j = nde_nax[F1]-ngst; j < nde_nax[F1]; ++j){

      int target_j = nde_nax[F1] - ngst -1 - (j - jRbnd);
      ntrack[j] = ntrack[target_j];
      nact[j] = nact[target_j];
      iRbnd[j] = iRbnd[target_j];
      iLbnd[j] = iLbnd[target_j];
      std::copy_n(&Ctot[target_j][0], ntrack[target_j],   &Ctot[j][0]);
      std::copy_n(&Itot[target_j][0], ntrack[target_j]-1, &Itot[j][0]);

      // updating ghost positions, ids and indexes
      for (int i = 0; i < ntrack[j]; ++i){
        int ind[] = {j,i};
        assignId(ind);
        Ctot[j][i].G.x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
        Ctot[j][i].computeAllGeom();
        if (i != ntrack[j]-1){ 
          Itot[j][i].x[F1] += (j-jRbnd+1)*Ctot[jRbnd-1][iLbnd[j]+1].G.dx[F1];
          Itot[j][i].computedA();
        }
      }
      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[UU2] *= -1;  
      }
    }
  }

  UNUSED(it);
  UNUSED(t);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell c = Ctot[j][i];                      // get current cell (track j, cell i)
  // double trac = c.S.prim[TR1];           // get tracer value
  double r   = c.G.x[r_];                   // get cell radial coordinate
  double dr  = c.G.dx[r_];                  // get cell radial spacing
  double dth = c.G.dx[t_];                  // get cell angular spacing
  
  Cell cOut = Ctot[j][iRbnd[j]];            // get outermost cell in track j
  double Rout = cOut.G.x[r_];               // get radial coordinate of outermost cell

  // determine which radius to use for AR calculations
  double r_AR;
  if (Rout < R2){                           // if furthest cell in track is within R2
    r_AR = Rout;
  }
  else if (r < R2){                         // if furthest cell in track is greater than R2 but cell is within R2
    r_AR = R2;
  }
  else{                                     // if cell is beyond R2
    r_AR = r;
  }

  // double AR = dr / (r * dth);               // calculate aspect ratio (only true AR if r_AR = r)
  double AR = dr / (r_AR * dth);            // calculate aspect ratio (only true AR if r_AR = r)
  
  // calculate target aspect ratio
  // double target_AR  = fmax(pow(R1/r, 0.2), 1.)*fin_AR;
  double target_AR  = fmax(pow(R1/r_AR, 0.2), 1.)*fin_AR;

  if (AR > split_AR * target_AR) {          // if cell is too long for its width
      return(split_);                       // split
  }
  if (AR < merge_AR * target_AR) {          // if cell is too short for its width
      return(merge_);                       // merge
  }
  return(skip_);
}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value
  
  UNUSED(*res);

}


void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(rho);
  UNUSED(p);
  UNUSED(uu);
  
  return;

}














