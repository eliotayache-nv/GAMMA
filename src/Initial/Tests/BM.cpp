#include "../../environment.h"
#include "../../grid.h"
#include "../../constants.h"

#ifndef BM_
#define BM_
#define LORENTZ_ 1
#endif

bool onedim = true;

// set shell and CBM parameters
static double n0      = 1.;           // cm-3:    CBM number density
static double rho0    = n0*mp_;       // g.cm-3:  comoving CBM mass density
static double eta     = 1.e-5;        //          eta = p/(rho*c^2)
static double th_simu = PI/32.;          // rad:     simulation angle
static double th_min  = 99.*PI/3200.;  //0. // rad:     minimum angle if avioding jet axis

// normalisation constants:
static double rhoNorm = rho0;                 // density normalised to CBM density
static double lNorm = c_;                     // distance normalised to c
static double vNorm = c_;                     // velocity normalised to c
static double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

// BM parameters
static double Etot = 1.e53;     // erg
static double n_ext = 1.e0;      // external medium number density
static double tinit = 4.e6;    // s,starting time (determines initial position of BW)
static double lfacstart = 1000;  // starting lfac behind the shock (determines initial position of BW)
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
  par->ncell[x_] = 100;              // number of cells in r direction
  par->ncell[y_] = 100;               // number of cells in theta direction
  if (onedim) par->ncell[y_] = 1;
  par->nmax      = 110;              // max number of cells in MV direction
  par->ngst      = 2;                 // number of ghost cells (?); probably don't change

  normalizeConstants(rhoNorm, vNorm, lNorm);

}

int Grid::initialGeometry(){                              

  // slighlty moving the grid to resolve the shock more precisely
  double rmin = rmin0;
  double rmax = rmax0;
  // double rmin = rmin0 + (rmax0-rmin0)/(2.*ncell[r_]);
  // double rmax = rmax0 + (rmax0-rmin0)/(2.*ncell[r_]);

  // grid defined in length units of light-seconds
  rmin /= lNorm;
  rmax /= lNorm;
  double dr = (rmax-rmin)/ncell[x_];

  // double logrmin  = log(rmin);
  // double logrmax  = log(rmax);
  // double dlogr    = (logrmax-logrmin);
  
  for (int j = 0; j < ncell[y_]; ++j){
    for (int i = 0; i < ncell[x_]; ++i){
      Cell *c = &Cinit[j][i];
      
      // double rlog   = (double) dlogr*(i+0.5)/ncell[x_] + logrmin;
      // double rlogL  = (double) dlogr*(i  )/ncell[x_] + logrmin;
      // double rlogR  = (double) dlogr*(i+1)/ncell[x_] + logrmin;
      // double r_ls   = exp(rlog);
      // double dr_ls  = exp(rlogR) - exp(rlogL);

      double r      = rmin + (i+0.5)*dr;
      double th     = (double) th_min + (th_simu - th_min)*(j+0.5)/ncell[y_];
      double dth    = (th_simu - th_min)/ncell[y_];

      c->G.x[x_]    = r;
      c->G.dx[x_]   = dr;
      c->G.x[y_]    = th;
      c->G.dx[y_]   = dth;

      c->computeAllGeom();
    }
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
  for (int j = 0; j < ncell[F1]; ++j){        // loop through cells along theta
    for (int i = 0; i < ncell[MV]; ++i){      // loop through cells along r
      Cell *c = &Cinit[j][i];

      double r = c->G.x[x_];                  // ls:  radial coordinate
      double dr = c->G.dx[x_];                  // ls:  radial coordinate
      double r_denorm = r*lNorm;
      double dr_denorm = dr*lNorm;

      if ( r_denorm > RShock){   // if in the shell tail
        double rho = n_ext*mp_*pow(r_denorm/Rscale, -k);
        c->S.prim[RHO]  = rho / rhoNorm;
        c->S.prim[VV1] = 0;
        c->S.prim[VV2] = 0;
        c->S.prim[PPP] = eta*rho*c_*c_/pNorm;
        c->S.prim[TR1] = 1.;
      }
      else{

        // computing an average over each cell
        int nbins = 1000;
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
        c->S.prim[VV2] = 0;
        c->S.prim[PPP] = p/pNorm;
        c->S.prim[TR1] = 2.;

        // c->S.prim2cons(c->G.x[r_]);
        // c->S.cons2prim(c->G.x[r_]);
      }
    }
  } 
  return 0;

}


void Grid::userKinematics(int it, double t){

  // setting lower and higher i boundary interface velocities
  // set by boundary velocity:
  // double vIn  = 0;
  double vOut = 1.05;

  int j = jLbnd+1;
  int i = iRbnd[j]-10;
  if (Ctot[j][i].S.prim[UU1] < 1.e-3 and it > 1000){
    vOut = 0.;  
  }
  
  for (int j = 0; j < nde_nax[F1]; ++j){
    for (int n = 0; n < ngst; ++n){
      int    iL = n;
      int    iR = ntrack[j]-2-n;
      // Itot[j][iL].v = vIn;
      Itot[j][iR].v = vOut;
    }
  }
}

void Cell::userSourceTerms(double dt){

  UNUSED(dt);

}

void Grid::userBoundaries(int it, double t){

  // for (int j = 0; j < nde_nax[F1]; ++j){
  //   for (int i = 0; i <= iLbnd[j]; ++i){
  //     Cell *c = &Ctot[j][i];
  //     double rho, u, p;
  //     double r = c->G.x[r_]*lNorm;
  //     calcBM(r, t, &rho, &u, &p);
  //     c->S.prim[RHO] = rho;
  //     c->S.prim[UU1] = u;
  //     c->S.prim[PPP] = p;
  //   }
  // }

  UNUSED(t);
  UNUSED(it);

}


int Grid::checkCellForRegrid(int j, int i){

  Cell c = Ctot[j][i];

  // double trac = c.S.prim[TR1];           // get tracer value
  double r   = c.G.x[r_];                   // get cell radial coordinate
  double dr  = c.G.dx[r_];                  // get cell radial spacing
  double dth = c.G.dx[t_];                  // get cell angular spacing
  double ar  = dr / (r*dth);                // calculate cell aspect ratio
  // double rOut = Ctot[j][iRbnd[j]-1].G.x[x_];
  // double ar  = dr / (rOut*dth);                // calculate cell aspect ratio

  // printf("%le\n", ar);

  int jtrk = jLbnd+1;
  double rmin = Ctot[jtrk][iLbnd[jtrk]+1].G.x[r_];
  double rmax = Ctot[jtrk][iRbnd[jtrk]-1].G.x[r_];
  double delta_r = rmax - rmin;
  double target_dr = delta_r / nact[jtrk];

  double split_DR   = 5.;                   // set upper bound as ratio of target_AR
  double merge_DR   = 0.3;                  // set upper bound as ratio of target_AR

  if (c.S.lfac() < 2.){
    split_DR = 20;
    merge_DR = 1.2;
  }

  if (dr > split_DR * target_dr) {          // if cell is too long for its width
    return(split_);                       // split
  }
  if (dr < merge_DR * target_dr) {          // if cell is too short for its width
    return(merge_);                       // merge
  }
  return(skip_);
}


void Cell::user_regridVal(double *res){
  // user function to find regrid victims
  // adapts the search to special target resolution requirements
  // depending on the tracer value

  if (S.lfac() < 2.){
    *res /= 4;
  }

}


void FluidState::cons2prim_user(double *rho, double *p, double *uu){

  UNUSED(uu);
  // double ratio = 1.e-5;
  // double p0   = eta*rho0*c_*c_;

  // *p = fmax(*p, ratio*p0/pNorm);
  // *rho = fmax(*rho, ratio*rho0/rhoNorm);

  return;

}














