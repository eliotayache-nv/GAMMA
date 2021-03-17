#ifndef CONSTANTS_H_
#define CONSTANTS_H_ 

#include "environment.h"

#define P_FLOOR_ 1.e-10         // pressure flooring

#define PI      M_PI            // Pi
#define mp_     1.6726219e-24   // Proton mass (g)
#define me_     9.1093835e-28   // electron mass (g)
#define qe_     4.80320425e-10  // electron charge (statC=cm3/2.g1/2.s-1) (rho01/2.l02.c)
#define c_      2.99792458e10   // speed of light (cm.s-1)
#define sigmaT_ 6.65345871e-25  // Thomson Cross section (cm2)
#define h_      6.62607004e-27  // Planck constant (cm2.g.s-1)

#define Msun_   1.98855e33      // Solar mass (g)

#define min_p_      1.e-10      // minimum value allowed for p
#define min_rho_    1.e-10      // minimum value allowed for rho

// Radiation related:
#define alpha_      1.29251816e-09  // alpha in gammaMax calculation (Van Eerten+2010)

// emissivity parameters:
#define p_          2.22    // slope of electron population (Kirk2000)
#define eps_e_      0.1     // constribution to electron acceleration
#define eps_B_      0.1     // contribution to magnetic field
#define zeta_       0.1    // fraction of accelerated electrons

// Normalized constants:
extern double Nmp_, Nme_, Nqe_, NsigmaT_, Nalpha_;

void normalizeConstants(double rhoNorm, double vNorm, double lNorm);

#endif
