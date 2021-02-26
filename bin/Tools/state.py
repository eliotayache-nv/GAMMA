# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2021-02-26 10:18:05
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-02-26 10:19:01

GAMMA_  = (5./3.)

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