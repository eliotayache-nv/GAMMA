#include "../../environment.h"
#include "../../grid.h"
#include "../../constants.h"
#include "../../simu.h"



// Structured jet params (Ryan, VE+ 2019)
static double log10_E0 = 52.96;
static double theta_c = 0.066;
static double theta_w = 0.47;
static double log10_n0 = -2.70;

// set shell and CBM parameters
static double n0      = pow(10, log10_n0);           // cm-3:    CBM number density
static double rho_ext = n0*mp_;       // g.cm-3:  comoving CBM mass density
static double eta     = 1.e-5;        //          eta = p/(rho*c^2)
static double th_simu = PI/2.;          // rad:     simulation angle
static double th_min  = 0.;//99.*PI/3200.; // rad:     minimum angle if avioding jet axis

// normalisation constants:
static double rhoNorm = rho_ext;                 // density normalised to CBM density
static double lNorm = c_;                     // distance normalised to c
static double vNorm = c_;                     // velocity normalised to c
static double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

// BM parameters
static double E0 = pow(10, log10_E0);     // erg
static double n_ext = n0;      // external medium number density
static double tstart = 3.3e7;    // s,starting time (determines initial position of BW) (4.36 for lfacf=100)
static double Rscale = 1.e17;   // cm
static double k = 0.;           // ext medium density profile

// intermediate quantities (Shock description)
static double lfacShock20 = E0 * (17.- 4.*k)/(8.*PI*n_ext*mp_*pow(tstart,3)*pow(c_,5));
  // Blandford&McKee(1976) eq. 69
static double RShock0 = c_*tstart*(1.-1./(1.+2.*(4.-k)*lfacShock20));
  // Blandford&McKee(1976) eq. 27

// grid size from shock position
static double rmin0 = RShock0*(1.-200./lfacShock20); // 1.5e6*lNorm;
static double rmax0 = RShock0*(1.+50./lfacShock20);

// ARM parameters
static double target_ar = 1.;     // target aspect ratio for lfac=1
static double split_AR   = 3.;    // set upper bound as ratio of target_AR
static double merge_AR   = 0.1;  // set upper bound as ratio of target_AR


static void calcBM(double E0, double n0,
                   double r, double t, double *rho, double *u, double *p, 
                   double *gmax, double *gmin){
  // returns normalised primitive variables for BM solution at position r (denormalised).

  // repeating from above (these are local)
  double lfacShock2 = E0 * (17.- 4.*k)/(8.*PI*n0*mp_*pow(t,3)*pow(c_,5));
    // Blandford&McKee(1976) eq. 69
  double RShock = c_*t*(1.-1./(1.+2.*(4.-k)*lfacShock2));
    // Blandford&McKee(1976) eq. 27
  double rhoa = n0*mp_*pow(RShock/Rscale, -k);
  double pa = rhoa*eta*c_*c_;
  FluidState Sa;
  Sa.prim[RHO] = rhoa;
  Sa.prim[PPP] = pa;
  double gma = Sa.gamma();  // Synge EOS
  double ha = c_*c_ + pa*gma/(gma-1.)/rhoa;
  double pf = 2./3.*lfacShock2*ha*rhoa;
  double lfacf = sqrt(fmax(1.,.5*lfacShock2));
  double Df = 2.*lfacShock2*rhoa;
    // Blandford&McKee(1976) eq. 8-10

  double chi = (1. + 2.*(4.-k)*lfacShock2) * (1. - r/(c_*t));
  *p = pf*pow(chi, -(17.-4.*k)/(12.-3.*k)) / pNorm;
  double lfac = sqrt(lfacf*lfacf/chi + 1);  // (+1) to ensure lfac>1
  double D = Df*pow(chi, -(7.-2.*k)/(4.-k));
    // Blandford&McKee(1976) eq. 28-30 / 65-67
  *rho = D/lfac / rhoNorm;
  *u = c_*sqrt(1.-1./(lfac*lfac))*lfac / vNorm;

  double dt0 = tstart / pow(chi, 0.25);
  double p0 = pf / pNorm;             // normalised because using EOS afterwards
  double rho0 = Df / lfacf / rhoNorm;          
  FluidState S0;
  S0.prim[RHO] = rho0;
  S0.prim[PPP] = p0;
  double gma0 = S0.gamma();
  double h0 = 1.+p0*gma0/(gma0-1.)/rho0;
  double eps0 = rho0*(h0-1.)/gma0;
  double eB0 = eps_B_ * eps0;
  double B0 = sqrt(8*PI*eB0);
  *gmax = 2.*19.*PI*Nme_*lfacf/(NsigmaT_*B0*B0*dt0) 
    * pow(chi, 25./24.)/(pow(chi, 19./12.)-1.); //gmax = +ifnty for chi=1
  *gmax = fmax(1.,*gmax);

  double psyn = p_;
  double ee0 = eps_e_ * eps0;
  double ne0 = zeta_ * rho0 / Nmp_;
  double lfac_av0 = ee0 / (ne0 * Nme_);
  double gmin0 = (psyn-2.) / (psyn-1.) * lfac_av0;
  *gmin = gmin0 / (pow(chi, 13./24.) + gmin0/ *gmax);
  *gmin = fmax(1.,*gmin);


}


void loadParams(s_par *par){

  par->tini      = tstart;             // initial time
  par->ncell[x_] = 9900;              // number of cells in r direction
  par->ncell[y_] = 300;               // number of cells in theta direction
  par->nmax      = 10000;              // max number of cells in MV direction
  par->ngst      = 2;                 // number of ghost cells (?); probably don't change

  normalizeConstants(rhoNorm, vNorm, lNorm);

}

int Grid::initialGeometry(){                              

  // slighlty moving the grid to resolve the shock more precisely
  // double dr0 = (rmax0-rmin0)/ncell[x_];
  // double shift = dr0 * (1 - (int) ((RShock - rmin0)/dr0)) * 0.99;
  // double rmin = rmin0 + shift;
  // double rmax = rmax0 + shift;
  double rmin = rmin0;
  double rmax = rmax0;

  // grid defined in length units of light-seconds
  rmin /= lNorm;
  rmax /= lNorm;

  // double logrmin  = log(rmin);
  // double logrmax  = log(rmax);
  // double dlogr    = (logrmax-logrmin);

  for (int j = 0; j < ncell[y_]; ++j){

    double dr_default = (rmax-rmin)/ncell[x_];
    double r_prec = rmin;
    double dr_prec = dr_default;

    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
            
      double dr     = dr_default;
      if (r_prec > RShock0/lNorm) { dr *= 10; } // increasing size of external cells
      double r      = r_prec + 0.5*dr_prec + 0.5*dr;
      double jrL    = (double)(j)/ncell[y_];
      double jrR    = (double)(j+1)/ncell[y_];
      double thL    = (double) th_min + (th_simu - th_min)*(0.3*jrL + 0.7*pow(jrL,3));
      double thR    = (double) th_min + (th_simu - th_min)*(0.3*jrR + 0.7*pow(jrR,3));
      double dth    = thR-thL;
      double th     = (thR+thL)/2.;

      c->G.x[x_]    = r;
      c->G.dx[x_]   = dr;
      c->G.x[y_]    = th;
      c->G.dx[y_]   = dth;

      c->computeAllGeom();

      r_prec = r;
      dr_prec = dr;
    }
  }

  return 0;
}

int Grid::initialValues(){

  // Careful, geometrical quantities are already normalised when using this function.
  
  if (GAMMA_ != 5./3.){
    printf("WARNING - Set GAMMA_ to 5./3.\n");
  }
  if (NUM_D != 2){
    printf("WARNING - set NUM_D to 2\n");
  }
  if (VI != 1.){
    printf("WARNING - Set VI to 1.\n");
  }
  if (SHOCK_DETECTION_ != ENABLED_){
    printf("WARNING - Set SHOCK_DETECTION_ to ENABLED_\n");
  }
  if (LOCAL_SYNCHROTRON_ != ENABLED_){
    printf("WARNING - Set LOCAL_SYNCHROTRON_ to ENABLED_\n");
  }
  if (eps_e_ != 0.038018){
    printf("Error! wrong initial value for eps_e_\n");
    exit(24);
  }
  if (eps_B_ != 0.000109){
    printf("Error! wrong initial value for eps_B_\n");
    exit(24);
  }
  if (zeta_ != 1.){
    printf("Error! wrong initial value for zeta_\n");
    exit(24);
  }


  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){        // loop through cells along theta
    for (int i = 0; i < ncell[MV]; ++i){      // loop through cells along r
      Cell *c = &Cinit[j][i];

      double r = c->G.x[x_];                  // ls:  radial coordinate
      double th = c->G.x[t_];                  // ls:  radial coordinate
      double r_denorm = r*lNorm;

      double theta2 = th*th;
      double theta_c2 = theta_c*theta_c;
      double E_track = E0*exp(-theta2/(2*theta_c2));

      double lfacShock2 = E_track * (17.- 4.*k)/(8.*PI*n_ext*mp_*pow(tstart,3)*pow(c_,5));
        // Blandford&McKee(1976) eq. 69
      double RShock = c_*tstart*(1.-1./(1.+2.*(4.-k)*lfacShock2)); // denormalised
        // Blandford&McKee(1976) eq. 27

      double rho, u, p, gmax, gmin;

      if ( r_denorm > RShock or th > theta_w){   // if in the shell tail
        rho = n_ext*mp_*pow(r_denorm/Rscale, -k);
        c->S.prim[RHO]  = rho / rhoNorm;
        c->S.prim[UU1] = 0;
        c->S.prim[UU2] = 0;
        c->S.prim[PPP] = eta*rho*c_*c_/pNorm;
        c->S.prim[TR1] = 1.;
        c->S.prim[GMX] = pow(rho, 1./3.);
        c->S.prim[GMN] = pow(rho, 1./3.);
        c->S.prim[PSN] = p_;
      }
      else{
        calcBM(E_track, n0, r_denorm, tstart, &rho, &u, &p, &gmax, &gmin);
        double v = u / sqrt(1.+u*u);
        c->S.prim[RHO] = rho;
        c->S.prim[VV1] = v;
        c->S.prim[VV2] = 0;
        c->S.prim[PPP] = p;
        c->S.prim[TR1] = 2.;
        c->S.prim[GMX] = pow(rho, 1./3.)/gmax;
        c->S.prim[GMN] = pow(rho, 1./3.)/gmin;
        c->S.prim[PSN] = p_;
      }

      // don't use cons2prim because doesn't have UU yet, but VV
      // c->S.prim2cons(c->G.x[r_]);
      // c->S.cons2prim(c->G.x[r_]);
    }
  } 

  if (worldrank == 0){
    printf("\n");
    printf("- 170817 SETUP -\n");
    printf("----------------\n");
    printf("\n");
    printf("Initial parameters:\n");
    printf("\n");
    printf("BM params:\n");
    printf("  E0      = %.1le (erg)\n", E0);
    printf("  n0      = %.1le (cm-3)\n", n0);
    printf("  tstart  = %.3le (s)\n", tstart);
    printf("\n");
    printf("  lfacS0  = %.0lf\n", sqrt(lfacShock20));
    printf("  RShock0 = %.4le\n", RShock0);
    printf("\n");
    printf("Jet params:\n");
    printf("  theta_c = %.2le (rad)\n", theta_c);
    printf("  theta_w = %.2le (rad)\n", theta_w);
    printf("  eta     = %.1le\n", eta);
    printf("\n");
    printf("Grid params:\n");
    printf("  nr     = %d \n  ntheta = %d\n", ncell[x_], ncell[y_]);
    printf("  radial size: %.4le (cm) -- %.4le (cm)\n", rmin0, rmax0);
    int nzones_in_BW = RShock0/lfacShock20 / ((rmax0 - rmin0)/ncell[r_]);
    printf("  blast-wave resolution: %d zones\n", nzones_in_BW);
    printf("\n");
    printf("AMR:\n");
    printf("  aspect ratio %lg -> %lg <- %lg\n", merge_AR*target_ar, 
                                                 target_ar, 
                                                 split_AR*target_ar);
    printf("\n");
    printf("----------------\n");
    printf("\n");
  }


  return 0;

}


void Grid::userKinematics(int it, double t){

  // setting lower and higher i boundary interface velocities
  // set by boundary velocity:
  // DEFAULT VALUES
  double vIn  = 0.;
  double vOut = 1.5;
  double vb = vOut;

  // Checking if shock is too far behind.
  if (worldrank == 0){
    int ja = jLbnd+1;
    int iout = iRbnd[ja]-1;
    int iin = iLbnd[ja]+1;
    double rout = Ctot[ja][iout].G.x[r_];
    double rin = Ctot[ja][iin].G.x[r_];
    double rlim = rin + 0.9 * (rout-rin);
    int ia = iin;
    double rcand = rin;
    Cell c = Ctot[ja][ia];
    double pcand = c.S.prim[PPP];
    while (rcand < rlim){
      ia++;
      c = Ctot[ja][ia];
      rcand = c.G.x[r_];
      pcand = c.S.prim[PPP];
    }
    if (pcand < 1.1*eta and t > 5.e6){ vb = 0; }
  }

  // applying updated bounary velocity to ll processes
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&vb, &vOut, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n < ngst; ++n){
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

  // REFLECTIVE PI/2 BOUNDARY
  if (worldrank == worldsize-1){
    for (int j = jRbnd; j < nde_nax[F1]; ++j){
      for (int i = 0; i <= ntrack[j]; ++i){
        Cell *c = &Ctot[j][i];
        c->S.prim[UU2] *= -1;
      }
    }
  }

  // BM BOUNDARY
  if (t<1.e7){
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i <= iLbnd[j]; ++i){
        Cell *c = &Ctot[j][i];
        double rho, u, p, gmax, gmin;
        double r = c->G.x[r_]*lNorm;
        double th = c->G.x[t_];

        double theta2 = th*th;
        double theta_c2 = theta_c*theta_c;
        double E_track = E0*exp(-theta2/(2*theta_c2));

        if (th < theta_w){
          calcBM(E_track, n0, r, t, &rho, &u, &p, &gmax, &gmin);
          c->S.prim[RHO] = rho;
          c->S.prim[UU1] = u;
          c->S.prim[UU2] = 0;
          c->S.prim[PPP] = p;
          c->S.prim[TR1] = 2.;
          c->S.prim[GMX] = pow(rho, 1./3.)/gmax;
          c->S.prim[GMN] = pow(rho, 1./3.)/gmin;
          c->S.prim[PSN] = p_;
        }
        else{
          c->S.prim[RHO]  = 1;
          c->S.prim[UU1] = 0;
          c->S.prim[UU2] = 0;
          c->S.prim[PPP] = eta;
          c->S.prim[TR1] = 1.;
          c->S.prim[GMX] = pow(1, 1./3.);
          c->S.prim[GMN] = pow(1, 1./3.);
        }
      }
    }
  }

  UNUSED(it);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell c = Ctot[j][i];

  // double trac = c.S.prim[TR1];           // get tracer value
  // double r   = c.G.x[r_];                   // get cell radial coordinate
  double dr  = c.G.dx[r_];                  // get cell radial spacing
  double dth = c.G.dx[t_];                  // get cell angular spacing
  // double ar_local  = dr / (r*dth);          // calculate cell aspect ratio
  // printf("%le\n", ar);
  double rOut = Ctot[j][iRbnd[j]-1].G.x[x_];
  double ar  = dr / (rOut*dth);             // calculate cell aspect ratio
                               
  double lfac = c.S.lfac();
  if (ar > split_AR * target_ar / lfac) { // if cell is too long for its width
    return(split_);                       // split
  }
  if (ar < merge_AR * target_ar / lfac) { // if cell is too short for its width
    return(merge_);                       // merge
  }
  return(skip_);
  
}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value
  
  double r = G.x[r_];
  double dr = G.dx[r_];
  double dth = G.dx[t_];
  double lfac = S.lfac();
  *res = dr / (r*dth) * lfac;
    // *lfac to give higher score to very small very fast cells

}


void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  // double ratio = 1.e-5;
  // double p0   = eta*rho0*c_*c_;

  // *p = fmax(*p, ratio*p0/pNorm);
  // *rho = fmax(*rho, ratio*rho0/rhoNorm);

  return;

}


void Simu::dataDump(){

  // if (it%1 == 0){ grid.printCols(it, t); }
  if (it%500 == 0){ grid.printCols(it, t); }

  // datadump in log time:
  int ndumps_per_decade = 1000;
  static double t_last_dump = tstart;
  double logdiff = log10(t)-log10(t_last_dump);
  if (logdiff > 1./ndumps_per_decade){
    t_last_dump = t; 
    grid.printCols(it, t); 
  }

}

void Simu::runInfo(){

  // if ((worldrank == 0) and (it%1 == 0)){ printf("it: %ld time: %le\n", it, t);}
  if ((worldrank == 0) and (it%100 == 0)){ printf("it: %ld time: %le\n", it, t);}

}

void Simu::evalEnd(){

  // if (it > 300){ stop = true; }
  if (t > 3.33e9){ stop = true; } // 3.33e8 BOXFIT simu

}












