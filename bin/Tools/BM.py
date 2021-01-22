# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2019-04-09 11:22:46
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-01-22 10:08:58

# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2019-03-26 16:53:26
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2019-04-09 11:02:34

# ---------------------------------------------------------------------------------------

import numpy as np

# ---------------------------------------------------------------------------------------

# Constants
#define PI      M_PI            // Pi
#define mp_     1.6726219e-24   // Proton mass (g)
#define me_     9.1093835e-28   // electron mass (g)
#define qe_     4.80320425e-10  // electron charge (statC=cm3/2.g1/2.s-1) (rho01/2.l02.c)
#define c_      2.99792458e10   // speed of light (cm.s-1)
#define sigmaT_ 6.65345871e-25  // Thomson Cross section (cm-2)
#define h_      6,62607004e-27  // Planck constant (cm2.g.s-1)

#define Msun_   1.98855e33      // Solar mass (g)

#define min_p_      1.e-10      // minimum value allowed for p
#define min_rho_    1.e-10      // minimum value allowed for rho

#define p_          2.2     // slope of electron population
#define eps_e_      0.1     // constribution to electron acceleration
#define eps_B_      0.1     // contribution to magnetic field

mp_ = 1.6726219e-24
c_  = 2.99792458e10
BM_p_rho_ratio_ = 1.e-3    # ratio between p and rho in ISM
GAMMA_ = (4./3.)


# ---------------------------------------------------------------------------------------


# TODO!
# Careful, c = 1 in all the calculations


class State(object):
    """docstring for State"""
    def __init__(self):
        super(State, self).__init__()
        
        self.rho = 0
        self.p = 0
        self.v = 0

        self.D = 0
        self.m = 0
        self.E = 0

        self.lfac = 0
        self.h    = 0

    def prim2aux(self):
        self.lfac = np.sqrt(1. / (1. - self.v**2 / c_**2))
        self.h = c_**2 + self.p * GAMMA_ / (GAMMA_ - 1) / self.rho


class BM(object):
    """docstring for BM"""

    def __init__(self, E, n_ext, t):
        super(BM, self).__init__()
        self.E          = E
        self.n_ext      = n_ext
        self.k          = 0
        self.t          = t

        self.rho_ext = self.n_ext * mp_
        self.lfacShock_sqrd = self.E * (17. - 4. * self.k) \
                                / (8. * np.pi * self.rho_ext * self.t**3 * c_**5)
            # Blandford&Mckee(1976) eq. 69
        self.lfacShock = np.sqrt(self.lfacShock_sqrd)

        print("lfacShock = %f, time = %f") %(self.lfacShock, self.t)

        self.RShock = self.t * (1. - 1. / (1. + 8. * self.lfacShock_sqrd)) * c_

        # Setting up front state
        self.Sf = State()    # Front
        self.Sa = State()    # Ahead of Front

        self.Sa.rho = self.rho_ext
        self.Sa.v = 0
        self.Sa.p = self.Sa.rho * BM_p_rho_ratio_ * c_**2

        self.Sa.prim2aux()

        self.Sf.p = 2./3. * self.lfacShock_sqrd * self.Sa.h * self.Sa.rho
        self.Sf.lfac = np.sqrt(max(1, 1./2. * self.lfacShock_sqrd))
        self.Sf.D = 2. * self.lfacShock_sqrd * self.Sa.rho
            # Blandford&Mckee(1976) eq. 8-10

        self.Sf.rho = self.Sf.D / self.Sf.lfac
        self.Sf.v = c_ * np.sqrt(1. - 1. / self.Sf.lfac)



    def local(self, r):

        S = State()

        if r > self.RShock:
            S.lfac = 1
            S.D = self.Sa.rho / (self.n_ext * mp_)
            S.p = self.Sa.p

        else:
            chi = (1. + 2. * (4.-self.k) * self.lfacShock_sqrd) * ( 1. - r / (c_ * self.t))
            S.p = self.Sf.p * chi**(-17. / 12.) 
            S.lfac = np.sqrt(self.Sf.lfac**2 / chi + 1)
            S.D = self.Sf.D * chi**(-7./4.)  / (self.n_ext * mp_)

        return(S)


def plotIsen1D(time, key, ax = None, **kwargs):

  U_ref = State(rho_ref,v_ref,p_ref)
  U_ref.prim2aux()
  U_ref.Jm = Jm(U_ref)
  U_ref.Jp = Jp(U_ref)

  # Evaluating initial conditions for reference:
  # These values are used to constrain the bisection used to derive future fluid states.
  x0, U0 = fillU0(U_ref)

  # Solving for later conditions:
  x = setup_x()

  # Computing exact solution:
  print("Computing exact solution:")
  U_arr = []
  for i in range(len(x)):
    if i%10 == 0:
      print("%i/%i" %(i, len(x)))
    U = findU(U_ref, x[i], time, U0)
    U_arr.append(U)

  if key == "rho":
    zth = np.array([U.rho for U in U_arr])
  if key == "vx":
    zth   = np.array([U.v for U in U_arr])
  if key == "p":
    zth   = np.array([U.p for U in U_arr])

  if ax is None:
    ax = plt.gca()

  ax.plot(x, zth, **kwargs)






















