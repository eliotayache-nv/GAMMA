#include "../../environment.h"
#include "../../grid.h"
#include "../../constants.h"
#include "../../simu.h"

#ifndef BM_
#define BM_
#define LORENTZ_ 0
#endif

bool onedim = true;
double dth = M_PI/2./1000.;


// set shell and CBM parameters
static double n0      = 1.;           // cm-3:    CBM number density
static double rho0    = n0*mp_;       // g.cm-3:  comoving CBM mass density
static double eta     = 1.e-5;        //          eta = p/(rho*c^2)

// normalisation constants:
static double rhoNorm = rho0;                 // density normalised to CBM density
static double lNorm = c_;                     // distance normalised to c
static double vNorm = c_;                     // velocity normalised to c
static double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

// BM parameters
static double Etot = 1.e53;     // erg
static double n_ext = 1.e0;      // external medium number density
static double tinit = 4.e6;    // s,starting time (determines initial position of BW)
static double lfacstart = 100;  // starting lfac behind the shock (determines initial position of BW)
static double Rscale = 1.e17;   // cm
static double k = 0.;           // ext medium density profile

// intermediate quantities (Shock description)
#if (LORENTZ_ == 1)
  static double lfacShock2 = 2*lfacstart*lfacstart;
  static double tstart = pow(Etot*(17.)/(8.*PI*n_ext* mp_*pow(c_,5.)*lfacShock2), 1./3.);
#else
  static double tstart = tinit;
  static double lfacShock2 = Etot * (17.- 4.*k)/(8.*PI*n_ext*mp_*pow(tstart,3)*pow(c_,5));
  // Blandford&McKee(1976) eq. 69}
#endif
static double RShock = c_*tstart*(1.-1./(1.+2.*(4.-k)*lfacShock2));
  // Blandford&McKee(1976) eq. 27
static double rhoa = n_ext*mp_*pow(RShock/Rscale, -k);
static double pa = rhoa*eta*c_*c_;
static double ha = c_*c_ + pa*GAMMA_/(GAMMA_-1.)/rhoa;
static double pf = 2./3.*lfacShock2*ha*rhoa;
static double lfacf = sqrt(fmax(1.,.5*lfacShock2));
static double Df = 2.*lfacShock2*rhoa;
  // Blandford&McKee(1976) eq. 8-10

// grid size from shock position
static double rmin0 = RShock*(1.-50./lfacShock2);
static double rmax0 = RShock*(1.+50./lfacShock2);

// ARM parameters
static double target_ar = 1.;     // target aspect ratio for lfac=1
static double split_AR   = 3.;    // set upper bound as ratio of target_AR
static double merge_AR   = 0.1;  // set upper bound as ratio of target_AR


static void calcBM(double r, double t, double *rho, double *u, double *p){
  // returns normalised primitive variables for BM solution at position r (denormalised).

  // repeating from above (these are local)
  double lfacShock2 = Etot * (17.- 4.*k)/(8.*PI*n_ext*mp_*pow(t,3)*pow(c_,5));
    // Blandford&McKee(1976) eq. 69
  double RShock = c_*t*(1.-1./(1.+2.*(4.-k)*lfacShock2));
    // Blandford&McKee(1976) eq. 27
  double rhoa = n_ext*mp_*pow(RShock/Rscale, -k);
  double pa = rhoa*eta*c_*c_;
  double ha = c_*c_ + pa*GAMMA_/(GAMMA_-1.)/rhoa;
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

}


void loadParams(s_par *par){

  par->tini      = tstart;             // initial time
  par->ncell[x_] = 9900;              // number of cells in r direction
  par->nmax      = 10000;              // max number of cells in MV direction
  par->ngst      = 2;                 // number of ghost cells (?); probably don't change

  normalizeConstants(rhoNorm, vNorm, lNorm);

}

int Grid::initialGeometry(){                              

  // slighlty moving the grid to resolve the shock more precisely
  double rmin = rmin0;
  double rmax = rmax0;
  rmin /= lNorm;
  rmax /= lNorm;
  double dr = (rmax-rmin)/ncell[x_];
  for (int i = 0; i < ncell[x_]; ++i){
    Cell *c = &Cinit[i];
    double r      = rmin + (i+0.5)*dr;
    c->G.x[x_]    = r;
    c->G.dx[x_]   = dr;
    c->computeAllGeom();
  }
  return 0;
}

int Grid::initialValues(){

  // Careful, geometrical quantities are already normalised when using this function.

  if (GAMMA_ != 4./3.){
    printf("WARNING - Set GAMMA_ to 4./3.\n");
  }
  if (VI != 1.){
    printf("WARNING - Set VI to 1.\n");
  }

  // initialise grid
  for (int i = 0; i < ncell[MV]; ++i){      // loop through cells along r
    Cell *c = &Cinit[i];

    double r = c->G.x[x_];                  // ls:  radial coordinate
    double dr = c->G.dx[x_];                  // ls:  radial coordinate
    double r_denorm = r*lNorm;
    double dr_denorm = dr*lNorm;

    if ( r_denorm > RShock){   // if in the shell tail
      double rho = n_ext*mp_*pow(r_denorm/Rscale, -k);
      c->S.prim[RHO]  = rho / rhoNorm;
      c->S.prim[VV1] = 0;
      c->S.prim[PPP] = eta*rho*c_*c_/pNorm;
      c->S.prim[TR1] = 1.;
    }
    else{

      // computing an average over each cell
      int nbins = 1;
      double ddx = dr_denorm / nbins;
      double xl = r_denorm - dr_denorm/2.;
      double xr = xl+ddx;
      double dV = 4./3.*PI*(pow(r_denorm+dr_denorm/2.,3) - pow(r_denorm-dr_denorm/2., 3));
      double rho=0,v=0,p=0;

      for (int ix = 0; ix < nbins; ++ix){
        double x = (xr+xl)/2.;
        double ddV = 4./3.*PI*(xr*xr*xr-xl*xl*xl);
        double chi = (1. + 2.*(4.-k)*lfacShock2) * (1. - x/(c_*tstart));
        double dp = pf*pow(chi, -(17.-4.*k)/(12.-3.*k));
        double dlfac = sqrt(lfacf*lfacf/chi + 1);  // (+1) to ensure lfac>1
        double dD = Df*pow(chi, -(7.-2.*k)/(4.-k));
          // Blandford&McKee(1976) eq. 28-30 / 65-67
        double drho = dD/dlfac;
        double dv = c_*sqrt(1.-1./(dlfac*dlfac));

        // printf("%le\n", v);

        rho += drho*ddV / dV;
        v += dv*ddV / dV;
        p += dp*ddV / dV;

        xl = xr;
        xr += ddx;
      }

      c->S.prim[RHO] = rho/rhoNorm;
      c->S.prim[VV1] = v/vNorm;
      c->S.prim[PPP] = p/pNorm;
      c->S.prim[TR1] = 2.;

    }
  } 
  return 0;

}


void Grid::userKinematics(int it, double t){

  // setting lower and higher i boundary interface velocities
  // set by boundary velocity:
  // double vIn  = 1.;
  // double vOut = 1.05;

  // int iout = iRbnd-10;
  // if (Ctot[iout].S.prim[UU1] < 1.e-3 and it > 1000){
  //   vOut = 0.;  
  // }
  // double xL = Ctot[iLbnd].G.x[r_];
  // double xR = Ctot[iRbnd].G.x[r_];
  // if (xL > 0.995*xR){
  //   vIn = 0.;
  // }
  
  // for (int n = 0; n < ngst; ++n){
  //   int    iL = n;
  //   int    iR = ntrack-2-n;
  //   Itot[iL].v = vIn;
  //   Itot[iR].v = vOut;
  // }

  // setting lower and higher i boundary interface velocities
  // set by boundary velocity:
  // DEFAULT VALUES
  double vIn  = 0.;
  double vOut = 1.5;
  double vb = vOut;

  // Checking if shock is too far behind.
  int iout = iRbnd-1;
  int iin = iLbnd+1;
  double rout = Ctot[iout].G.x[r_];
  double rin = Ctot[iin].G.x[r_];
  double rlim = rin + 0.9 * (rout-rin);
  int ia = iin;
  double rcand = rin;
  Cell c = Ctot[ia];
  double pcand = c.S.prim[PPP];
  while (rcand < rlim){
    ia++;
    c = Ctot[ia];
    rcand = c.G.x[r_];
    pcand = c.S.prim[PPP];
  }
  if (pcand < 1.1*eta and t > 5.e6){ vb = 0; }

  vb = vOut;
  
  for (int n = 0; n < ngst; ++n){
    int    iL = n;
    int    iR = ntrack-2-n;
    Itot[iL].v = vIn;
    Itot[iR].v = vOut;
  }

}

void Cell::userSourceTerms(double dt){

  UNUSED(dt);

}

void Grid::userBoundaries(int it, double t){

  if (t<1.e7){
    for (int i = 0; i <= iLbnd; ++i){
      Cell *c = &Ctot[i];
      double rho, u, p;
      double r = c->G.x[r_]*lNorm;

      calcBM(r, t, &rho, &u, &p);
      c->S.prim[RHO] = rho;
      c->S.prim[UU1] = u;
      c->S.prim[UU2] = 0;
      c->S.prim[PPP] = p;
      c->S.prim[TR1] = 2.;
    }
  }

  UNUSED(t);
  UNUSED(it);

}


int Grid::checkCellForRegrid(int j, int i){

  UNUSED(j);
  Cell c = Ctot[i];

  // double trac = c.S.prim[TR1];           // get tracer value
  // double r   = c.G.x[r_];                   // get cell radial coordinate
  double dr  = c.G.dx[r_];                  // get cell radial spacing
  // double ar_local  = dr / (r*dth);          // calculate cell aspect ratio
  // printf("%le\n", ar);
  double rOut = Ctot[iRbnd-1].G.x[x_];
  double ar  = dr / (rOut*dth);             // calculate cell aspect ratio

  // We will enforce stronger refinement closer to the shock front in the radial direction
  // Looking for the shock in the j track.
  double iS = iLbnd+1;
  for (int ii = iRbnd-1; ii > iLbnd; --ii){
    Cell ci = Ctot[ii];
    double p = ci.S.prim[PPP];
    iS = ii;
    if (p > 1.5*eta){ // The shock is showing on the left of cell ii
      // let's return the shock index in the track
      break;
    }
  }

  // Let's now check if we are close to the shock position
  double dist = i-iS;
  if (0 < dist and dist < 10){ ar *= 10; }  
                                             // -3 because we want the whole shock to be
                                             // resolved. this is hard-wired for now but
                                             // could be changed
                               
  double lfac = c.S.lfac();
  if (ar > split_AR * target_ar / pow(lfac, 3./2.)) { // if cell is too long for its width
    return(split_);                       // split
  }
  if (ar < merge_AR * target_ar / pow(lfac, 3./2.)) { // if cell is too short for its width
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
  double lfac = S.lfac();
  *res = dr / (r*dth) * pow(lfac, 3./2.);

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
  if (t > 3.33e8){ stop = true; } // 3.33e8 BOXFIT simu

}














