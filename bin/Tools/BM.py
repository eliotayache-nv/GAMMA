# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2019-04-09 11:22:46
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-01-27 17:59:02

# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2019-03-26 16:53:26
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2019-04-09 11:02:34

# ---------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------------

c_      = 2.99792458e10   # speed of light (cm.s-1)
mp_     = 1.6726219e-24   # Proton mass (g)
me_     = 9.1093835e-28   # electron mass (g)
sigmaT_ = 6.65345871e-25  # Thomson Cross section (cm2)
GAMMA_  = (4./3.)

BM_p_rho_ratio_ = 1.e-5    # ratio between p and rho in ISM
E0 = 1.e53
n0 = 1.e0
rho0 = n0*mp_

# Emissivity parameters
p = 2.3
eps_B_ = 0.1
eps_e_ = 0.1
zeta_ = 0.1

# Normalisation
lNorm = c_
vNorm = c_
rhoNorm = rho0
pNorm = rhoNorm * vNorm * vNorm


# ---------------------------------------------------------------------------------------
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

    self.gmax = 0

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
      S.D = self.Sa.rho
      S.p = self.Sa.p
      S.rho = S.D/S.lfac
      S.gmax = 1.

    else:
      chi = (1. + 2. * (4.-self.k) * self.lfacShock_sqrd) * ( 1. - r / (c_ * self.t))
      S.p = self.Sf.p * chi**(-17. / 12.) 
      S.lfac = np.sqrt(self.Sf.lfac**2 / chi + 1)
      S.D = self.Sf.D * chi**(-7./4.) 
      S.rho = S.D/S.lfac

      t0 = self.t/chi**(1./4.)
      p0 = self.Sf.p
      rho0 = self.Sf.rho
      h0 = c_**2 + p0*GAMMA_/(GAMMA_-1.)/rho0
      eps0 = rho0*(h0 - c_**2) / GAMMA_
      eB0 = eps_B_ *eps0
      B0 = np.sqrt(8.*np.pi*eB0)
      S.gmax = 2.*19.*np.pi*me_*c_*self.Sf.lfac / (sigmaT_*B0**2*t0) \
        * chi**(25./24.) / (chi**(19./12.)-1.)
      print(S.gmax)

    return(S)


def plotBM1D(data, key, ax = None, **kwargs):

  time = data["t"][0]
  x = data["x"]*lNorm
  y = np.zeros(x.shape[0])

  BW = BM(E0, n0, time)

  for i in range(x.shape[0]):
    if key == "rho":
      y[i] = BW.local(x[i]).rho/rhoNorm
    if key == "lfac":
      y[i] = BW.local(x[i]).lfac
    if key == "p":
      y[i] = BW.local(x[i]).p/pNorm
    if key == "gmax":
      y[i] = BW.local(x[i]).gmax

  if ax is None:
    ax = plt.gca()

  ax.plot(x/lNorm, y, **kwargs)





















