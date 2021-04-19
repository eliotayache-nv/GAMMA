#include "../environment.h"
#include "../grid.h"
#include "../constants.h"

bool onedim = true;

// set shell and CBM parameters
static double n0      = 1.;           // cm-3:    CBM number density
static double rho0    = n0*mp_;       // g.cm-3:  comoving CBM mass density
static double eta     = 1.e-5;        //          eta = p/(rho*c^2)
static double th_simu = PI/32.;          // rad:     simulation angle
static double th_min  = 99.*PI/3200.;  //0. // rad:     minimum angle if avioding jet axis
static double t_ini   = 0.;         // s:     start time of dynamic simulation

// set simulation boundary positions
static double rmin    = 1.;           // light-seconds, box lower boundary at t_ini
static double rmax    = 500.;         // light-seconds, box upper boundary at t_ini

// normalisation constants:
static double rhoNorm = rho0;                 // density normalised to CBM density
static double lNorm = c_;                     // distance normalised to c
static double vNorm = c_;                     // velocity normalised to c
static double pNorm = rhoNorm*vNorm*vNorm;    // pressure normalised to rho_CMB/c^2

// Fireball parameters
static double R0   = 50.;      // light-seconds
static double Etot = 1.e52;     // erg
static double M0   = 3.6e29;    // g


void loadParams(s_par *par){

  par->tini      = t_ini;             // initial time
  par->ncell[x_] = 200;              // number of cells in r direction
  par->ncell[y_] = 100;               // number of cells in theta direction
  if (onedim) par->ncell[y_] = 1;
  par->nmax      = 1000;              // max number of cells in MV direction
  par->ngst      = 2;                 // number of ghost cells (?); probably don't change

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
      
      double th     = (double) th_min + (th_simu - th_min)*(j+0.5)/ncell[y_];
      double dth    = (th_simu - th_min)/ncell[y_];

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

  double c2   = c_*c_;
  double p0   = eta*rho0*c2;                  // CBM pressure
  
  if (GAMMA_ != 4./3.){
    printf("WARNING - Set GAMMA_ to 4./3.\n");
  }
  if (VI != 1.){
    printf("WARNING - Set VI to 1.\n");
  }

  // declare variables for initial setup integrals
  // NOTE: These are isotropic values. For the actual injection energy and mass, multiply by jet_solid_angle/4pi
  double Eiso = 0.;                           // initialise for energy calculation (integral)
  double Miso = 0.;                           // initialise for mass calculation (integral)
  double f0   = 0.;
  double f1   = 0.;
  double m0   = 0.;
  double m1   = 0.;
  double r0   = 0.;
  double r1   = 0.;

  double V0 = 4./3.*M_PI*R0*R0*R0*c_*c_*c_; 
  double D0 = M0/V0;
  double E0 = Etot/V0;
  double tau0 = E0 - D0*c2;
  printf("%le %le\n", E0, tau0);

  // initialise grid
  for (int j = 0; j < ncell[F1]; ++j){        // loop through cells along theta
    for (int i = 0; i < ncell[MV]; ++i){      // loop through cells along r
      Cell *c = &Cinit[j][i];

      double r = c->G.x[x_];                  // ls:  radial coordinate
      double th = c->G.x[y_];                 // rad: theta angular coordinate
              
      if ( r <= R0){   // if in the shell tail          
        c->S.cons[DEN]  = D0/rhoNorm;
        c->S.cons[SS1]  = 0.;
        c->S.cons[SS2]  = 0.;
        c->S.cons[TAU]  = tau0/(rhoNorm * c2);
        c->S.prim[TR1] = 1.;
        c->S.cons2prim(r); // r is already normalised
        // printf("%le %le\n", c->S.prim[PPP], c->S.prim[RHO]);
      }

      else{
        c->S.prim[RHO] = rho0/rhoNorm;
        c->S.prim[VV1] = 0;
        c->S.prim[VV2] = 0;
        c->S.prim[PPP] = p0/pNorm;
        c->S.prim[TR1] = 2.;
      }
      
      
      // integrate energy; only do integral along 1 track, otherwise this adds up all tracks in the integral giving a factor of Ntheta higher
      if (j == 0){                                          // if in first track
        double dens = c->S.prim[RHO]*rhoNorm;
        double vel  = c->S.prim[VV1]*vNorm;
        double rad  = r*lNorm;
        double lfac = 1./sqrt(1.-((vel*vel)/(c_*c_)));
        double Edot = dens*4*PI*(rad*rad)*vel*(lfac*lfac)*(c_*c_)*(1.+(eta*((GAMMA_/(GAMMA_-1.)-(1./(lfac*lfac))))));
      
        if(i == 0){
          r1 = rad;
          m1 = lfac*dens*4*PI*rad*rad;
          if (Edot == 0){
            f1 = 0;
          }
          else{
            f1 = Edot/vel;
          }
        }
        else{
          r0 = r1;
          f0 = f1;
          m0 = m1;
          r1 = rad;
          m1 = lfac*dens*4*PI*rad*rad;
          if (Edot == 0){
            f1 = 0;
          }
          else{
            f1 = Edot/vel;
          }
          double dr = r1-r0;
          Eiso += 0.5*(f0 + f1)*dr;
          Miso += 0.5*(m0 + m1)*dr;
        }
      }
    }
  }
  printf("Eiso = %e\tMiso = %e\n",Eiso,Miso);
  return 0;
}


void Grid::userKinematics(int it, double t){

  // setting lower and higher i boundary interface velocities

  // set by boundary velocity:
  double vIn  = 0.75;
  double vOut = 1.05;


  // set by Lorentz factor:
  // double vIn_lfac   = 100; // inner boundary
  // double vOut_lfac  = 200; // outer boundary
  // double vIn        = sqrt(1.-(1./(vIn_lfac*vIn_lfac)))*c_/vNorm;
  // double vOut       = sqrt(1.-(1./(vOut_lfac*vOut_lfac)))*c_/vNorm;     
  
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

  Cell c = Ctot[j][i];

  // double trac = c.S.prim[TR1];           // get tracer value
  double r   = c.G.x[r_];                   // get cell radial coordinate
  double dr  = c.G.dx[r_];                  // get cell radial spacing
  double dth = c.G.dx[t_];                  // get cell angular spacing
  double ar  = dr / (r*dth);                // calculate cell aspect ratio
  // double rOut = Ctot[j][iRbnd[j]-1].G.x[x_];
  // double ar  = dr / (rOut*dth);                // calculate cell aspect ratio

  
  // aspect ratio based regridding:
  // calculate target aspect ratio at radius r:
  // double sim_time   = 3.e7;                 // set expected simulation end time
  // double start_AR   = 10.;                  // set initial target aspect ratio
  // double end_AR     = 3.;                   // set target aspect ratio at expected end time
  // double target_AR  = start_AR - ((start_AR-end_AR)/sim_time)*r;
  // the above does not work - gives low resolution between start and end points
  // empirically found the following to work:

  // double y0         = 0.9863;
  // double yM         = 70.23;
  // double k          = 1.660;
  // double target_AR  = (y0 + (yM-y0) * exp(-k*log10(r)));
  double target_AR  = 3;

  // new
  // double y0         = 1.0146;
  // double yM         = 3495.2;
  // double k          = 2.790;
  // double target_AR  = y0 + (yM-y0) * exp(-k*log10(r));

  //double target_AR  = 1.;

  double split_AR   = 5.;                   // set upper bound as ratio of target_AR
  double merge_AR   = 0.2;                  // set upper bound as ratio of target_AR

  if (ar > split_AR * target_AR) {          // if cell is too long for its width
      return(split_);                       // split
  }
  if (ar < merge_AR * target_AR) {          // if cell is too short for its width
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

  UNUSED(uu);
  // double ratio = 1.e-5;
  // double p0   = eta*rho0*c_*c_;

  // *p = fmax(*p, ratio*p0/pNorm);
  // *rho = fmax(*rho, ratio*rho0/rhoNorm);

  return;

}














