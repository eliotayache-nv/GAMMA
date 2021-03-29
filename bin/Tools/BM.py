# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2019-04-09 11:22:46
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-03-28 23:00:55

# ---------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import glob

from gamma_io import *

# ---------------------------------------------------------------------------------------

c_      = 2.99792458e10   # speed of light (cm.s-1)
mp_     = 1.6726219e-24   # Proton mass (g)
me_     = 9.1093835e-28   # electron mass (g)
sigmaT_ = 6.65345871e-25  # Thomson Cross section (cm2)
GAMMA_  = (5./3.)

BM_p_rho_ratio_ = 1.e-5    # ratio between p and rho in ISM
E0 = 1.e53
n0 = 1.e0
rho0 = n0*mp_

# Emissivity parameters
p_ = 2.22
eps_B_ = 0.1
eps_e_ = 0.1
zeta_ = 0.1

# Normalisation
lNorm = c_
vNorm = c_
rhoNorm = rho0
pNorm = rhoNorm * vNorm * vNorm

Nme_ = me_/(rhoNorm * lNorm**3)
Nmp_ = mp_/(rhoNorm * lNorm**3)
NsigmaT_ = sigmaT_ / lNorm**2


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
    self.gmin = 0

  def gamma(self):
    rho = self.rho
    p = self.p
    a = p/rho / (GAMMA_-1)
    e_ratio = a + np.sqrt(a**2+1.)
    gamma_eff = GAMMA_ - (GAMMA_-1.)/2. *(1.-1./(e_ratio**2))
    return gamma_eff
 
  def prim2aux(self):
    self.lfac = np.sqrt(1. / (1. - self.v**2 / c_**2))
    gma = self.gamma()
    self.h = c_**2 + self.p * gma / (gma - 1) / self.rho


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
      p0 = self.Sf.p / pNorm
      rho0 = self.Sf.rho /rhoNorm
      gma0 = self.Sf.gamma()
      h0 = 1. + p0*gma0/(gma0-1.)/rho0
      eps0 = rho0*(h0 - 1.) / gma0
      eB0 = eps_B_ * eps0
      B0 = np.sqrt(8.*np.pi*eB0)
      S.gmax = 2.*19.*np.pi*Nme_*self.Sf.lfac / (NsigmaT_*B0**2*t0) \
        * chi**(25./24.) / (chi**(19./12.)-1.) # gmax = +ifnty for chi=1

      psyn = p_
      ee0 = eps_e_ * eps0
      ne0 = zeta_ * rho0 / Nmp_
      lfac_av0 = ee0 / (ne0 * Nme_)
      gmin0 = (psyn-2.) / (psyn-1.) * lfac_av0
      S.gmin = gmin0 / (chi**(13./24.) + gmin0/ S.gmax)

    return(S)


def plotBM1D(data, key, jtrack=None, x_norm=None, ax = None, **kwargs):

  time = data["t"][0]

  if jtrack is not(None):
    x = pivot(data, "x")[jtrack,:]
    x = x[~np.isnan(x)]*lNorm
  else:
    x = np.copy(data["x"]*lNorm)
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
    if key == "gmin":
      y[i] = BW.local(x[i]).gmin

  if ax is None:
    ax = plt.gca()

  x_plot = np.copy(x)/lNorm
  if x_norm is not(None):
    x_plot /= x_norm

  ax.plot(x_plot, y, **kwargs)


def AsymptoteBM(dir="Last"):

  lfacMaxArr = []
  lfacThArr = []
  lfacSArr = []
  tArr = []
  rArr = []
  rThArr = []

  for filename in glob.glob("../../results/%s/*.out" %dir):
    print(filename.split("results/")[1])
    data = readData(filename.split("results/")[1])

    t = data["t"][0]
    r = data["x"].to_numpy()
    v = data["vx"].to_numpy()
    lfac = 1./np.sqrt(1 - v**2)
    lfacMax = np.max(lfac)
    imax = np.argmax(lfac)
    rS = r[imax]

    BW = BM(E0, n0, t)
    lfacf = BW.Sf.lfac
    lfacS = BW.lfacShock
    r_th = BW.RShock

    tArr.append(t)
    rArr.append(rS)
    rThArr.append(r_th)
    lfacMaxArr.append(lfacMax)
    lfacThArr.append(lfacf)
    lfacSArr.append(lfacS)

  lfac = np.array(lfacMaxArr)
  lfac_th = np.array(lfacThArr)
  lfacS_th = np.array(lfacSArr)
  t = np.array(tArr)   
  rS = np.array(rArr)
  rS_th = np.array(rThArr)

  tinds = t.argsort()
  tsort = t[tinds]
  rS = rS[tinds]
  vS = (rS[1:]-rS[:-1])/(tsort[1:]-tsort[:-1])
  lfacS = 1./np.sqrt(1.-vS**2)
  lfacS_th = lfacS_th[tinds]

  plt.plot(tsort, lfac[tinds])
  plt.plot(tsort, lfac_th[tinds])

  return(t, tsort, lfac, lfac_th, lfacS, lfacS_th, rS, rS_th)


def BM_precision(data, key):

  time = data["t"][0]
  x = data["x"]
  dx = data["dx"]
  y = data[key]
  y_th = np.zeros(x.shape[0])

  BW = BM(E0, n0, time)

  for i in range(x.shape[0]):
    if key == "rho":
      y_th[i] = BW.local(x[i]*lNorm).rho/rhoNorm
    if key == "lfac":
      y_th[i] = BW.local(x[i]*lNorm).lfac
    if key == "p":
      y_th[i] = BW.local(x[i]*lNorm).p/pNorm
    if key == "gmax":
      y_th[i] = BW.local(x[i]*lNorm).gmax

  L1 = computeL1(y, y_th, dx)
  return(L1)


def BM_convergence(key, dir="Last"):

  L1 = []
  t = []

  for filename in glob.glob("../../results/%s/*.out" %dir):
    print(filename.split("results/")[1])
    data = readData(filename.split("results/")[1])
    t.append(data["t"][0])
    L1.append(BM_precision(data, key))

  t = np.array(t)
  L1 = np.array(L1)

  plt.scatter(t, L1)

























