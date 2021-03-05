# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2021-02-26 10:15:40
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-03-04 16:21:33

from state import *

c_      = 2.99792458e10   # speed of light (cm.s-1)
mp_     = 1.6726219e-24   # Proton mass (g)
me_     = 9.1093835e-28   # electron mass (g)
sigmaT_ = 6.65345871e-25  # Thomson Cross section (cm2)
qe_     = 4.80320425e-10  # electron charge (statC=cm3/2.g1/2.s-1) (rho01/2.l02.c)

# Emissivity parameters
p_ = 2.3
eps_B_ = 0.1
eps_e_ = 0.1
zeta_ = 0.1

# Normalisation
n0 = 1.e0
rho0 = n0*mp_

lNorm = c_
vNorm = c_
rhoNorm = rho0
pNorm = rhoNorm * vNorm * vNorm

Nme_ = me_/(rhoNorm * lNorm**3)
NsigmaT_ = sigmaT_ / lNorm**2
Nqe_ = qe_ / (rhoNorm**0.5 * lNorm**2 * vNorm)


def nuPrime(nu, lfac, beta, mu):
  return(nu*lfac*(1.-beta*mu))

def nuSyn(gammae, B):
  return(3*gammae**2*Nqe_*B/(16.*Nme_))


def SynchPower(S, theta, phi, thetaobs, nuobs):

  rho = S.rho
  p = S.p
  gma = S.gamma()
  h = 1. + p*gma/(gma-1.)/rho
  eps = rho*(h - 1.) / gma
  eB = eps_B_ * eps
  B = np.sqrt(8.*np.pi*eB)

  gmax = S.gmax
  numax = nuSyn(gmax, B)
  gmin = S.gmin
  numin = nuSyn(gmin, B)

  vr   = S.v[0]
  vth  = S.v[1]
  beta = S.vel()
  lfac = S.lfac()
  thetav = theta + np.arctan2(vth,vr)  # angle between velocity and jet axis (phi stays same as cell)
  thetaeff = np.arccos(np.sin(thetav)*np.cos(phi)*np.sin(thetaobs) 
                       + np.cos(thetav)*np.cos(thetaobs))
  mu = np.cos(thetaeff)
  nu = nuPrime(nuobs, lfac, beta, mu)

  n = rho / Nme_;
  Pmax = 4.*(p_-1)/(3*p_-1) * n * NsigmaT_ \
    * 4./3. * B / (6.*np.pi) * 16.* Nme_ / (3.*Nqe_)

  if nu < numin:
    P = Pmax * (nu/numin)**(1./3.)
  if numin < nu < numax:
    P = Pmax * (nu/numin)**(-(p_-1.)/2.)
  if numax < nu:
    P = 0

  return(P)


def PSynFromPandas(rho, p, vr, vth, gmax, gmin, theta, phi, thetaobs, nuobs):
  S = State()  
  S.rho = rho
  S.p = p
  S.v[0] = vr
  S.v[1] = vth
  S.gmax = gmax
  S.gmin = gmin
  return(SynchPower(S, theta, phi, thetaobs, nuobs))


def addSynchPower(data, phi, thetaobs, nuobs):
  PSyn = [PSynFromPandas(rho, p, vr, vth, gmax, gmin, theta, phi, thetaobs, nuobs) 
    for rho, p, vr, vth, gmax, gmin, theta 
    in zip(data["rho"], data["p"], data["vx"], data["vy"], data["gmax"], data["gmin"], 
           data["y"])]
  return(np.array(PSyn))

