# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2021-02-26 10:18:05
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-03-01 20:21:23

import numpy as np

GAMMA_  = (5./3.)

class State(object):
  """docstring for State"""
  def __init__(self):
    super(State, self).__init__()
    
    self.rho = 0
    self.p = 0
    self.v = np.zeros(2)

    self.D = 0
    self.m = 0
    self.tau = 0

    self.gmax = 0
    self.gmin = 0

  def vel(self):
    return(np.sqrt(np.sum(self.v**2)))

  def lfac(self):
    return(1./np.sqrt(1.-self.vel()**2))

  # Synge type EOS
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