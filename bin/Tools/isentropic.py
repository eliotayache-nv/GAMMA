# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2019-12-13 11:34:21
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-12-07 17:29:58


"""
This program computes the exact solution of an isentropic wave evolution using the method
of characteristics and outputs the result to a file.
"""


# ----------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
import os

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rc('xtick', labelsize=18) 
# plt.rc('ytick', labelsize=18) 
# plt.rcParams['savefig.dpi'] = 300


# ----------------------------------------------------------------------------------------
# Constants
GAMMA_ = 5./3.
c_     = 1.


# ----------------------------------------------------------------------------------------
# Wave parameters
rho_ref = 1.
v_ref   = 0.
p_ref   = 1.e2
alpha   = 1.
L       = 0.3

# Domain parameters
xmin    = -0.35
xmax    = 1
npts    = 500

# Simulation parameters
tstart  = 0.
teval   = 1.5
x_prec  = 1.e-14    # ideally, 1.e-15 (Machine prec)


# ----------------------------------------------------------------------------------------
def parseArguments():
  """Parses the optional arguments to the program"""

  # Create argument parser 
  parser = argparse.ArgumentParser()

  parser.add_argument("-d","--datadir", 
    help="directory containing phys.out", type=str, default="../../results/Last")
  parser.add_argument("-f","--filename", 
    help="File name for comparison", type=str, default="phys0000000000.h5")
  parser.add_argument("-i","--initial-state", 
    help="Overplot the initial state",
    action='store_true')
  parser.add_argument("-p","--plot", 
    help="Do you want to plot the results?",
    action='store_true')
  parser.add_argument("-c","--cell-size", 
    help="Do you want to plot the local cell size as well?",
    action='store_true')

   # Parse arguments
  args = parser.parse_args()
  return args


# ----------------------------------------------------------------------------------------
def openfile(filename):
  """Open HDF5 files and check for fails"""
  try:
    file = h5py.File(filename, 'r')
    status = 0
    return(file, status)
  except IOError:
    print("This file doesn't exist. Please choose another file.")
    status = 1
    return(None, status)


# ----------------------------------------------------------------------------------------
class State(object):
  """Fluid state (primitive, conserved and auxilliary variables"""
  def __init__(self, rho=0, v=0, p=0):
    super(State, self).__init__()

    self.rho = rho
    self.v = v
    self.p = p

    self.D = 0.
    self.m = 0.
    self.E = 0.

    self.lfac = 0.
    self.h    = 0.
    self.cs   = 0.

    self.Jp   = 0.
    self.Jm   = 0.

  def prim2aux(self):
    """Computes lfac, h, cs, lm and lp from rho, v and p"""
    self.lfac = np.sqrt(1. / (1. - self.v**2 / c_**2))
    self.h    = c_**2 + self.p * GAMMA_ / (GAMMA_ - 1.) / self.rho
    self.cs   = np.sqrt((GAMMA_ * self.p) / (self.rho*self.h))
    self.lm   = (self.v - self.cs) / (1. - self.v*self.cs)   # lambda_minus
    self.lp   = (self.v + self.cs) / (1. + self.v*self.cs)   # lambda_plus


# ----------------------------------------------------------------------------------------
def rho0(x):
  """Initial mass density profile"""
  if abs(x) < L:
    f = ((x/L)**2 - 1.)**4
  else:
    f = 0
  return(rho_ref * (1 + alpha * f))


# ----------------------------------------------------------------------------------------
def Jp(U):
  """Returns the Riemann invariants JM corresponding to a given state vector"""
  Jp = 1./2. * np.log((1+U.v)/(1-U.v)) \
    + 1./np.sqrt(GAMMA_-1.) * np.log((np.sqrt(GAMMA_-1)+U.cs) / (np.sqrt(GAMMA_-1)-U.cs))
  return(Jp)

def Jm(U):
  """Returns the Riemann invariants JM corresponding to a given state vector"""
  Jm = 1./2. * np.log((1+U.v)/(1-U.v)) \
    - 1./np.sqrt(GAMMA_-1.) * np.log((np.sqrt(GAMMA_-1)+U.cs) / (np.sqrt(GAMMA_-1)-U.cs))
  return(Jm)


# ----------------------------------------------------------------------------------------
def computeU0(U_ref, x):
  """Return a initial state vcector from x and U_ref"""
  U = State(rho0(x))

  if abs(x) > L:
    U = U_ref
  else:
    U.Jm = U_ref.Jm
    U.p = U_ref.p * (U.rho / U_ref.rho)**GAMMA_
    U.prim2aux()
      # Careful! lfac is wrong there!
    temp = np.exp(2.*U.Jm) \
      * ((np.sqrt(GAMMA_-1.)+U.cs) / (np.sqrt(GAMMA_-1.)-U.cs))**(2./np.sqrt(GAMMA_-1.))
    U.v = (temp-1)/(temp+1)
    U.prim2aux()
      # lfac is set right here
    U.Jp = Jp(U)

  return(U)


# ----------------------------------------------------------------------------------------
def fillU0(U_ref):
  """Sets all fluid values to their initial state from x0 and rho0"""

  x0 = np.linspace(xmin,xmax,npts)

  U0 =[]
  for i in range(npts):
    U = computeU0(U_ref, x0[i])
    U0.append(U)
 
  return x0, U0


# ----------------------------------------------------------------------------------------
def setup_x(compare=False):
  """Sets the array of x at teval"""
  x = np.linspace(xmin, xmax, npts)
  return x


# ----------------------------------------------------------------------------------------
def findU(U_ref, x, teval, U0_arr):
  """
  Finds the value of U for a value of x, teval and U_ref
  """

  # We target the value of Jp, using the initial values of lp.
  # Jm and s and constant across all of spacetime.

  # Building initial bisection domain
  lmin = np.min(np.array([U.lp for U in U0_arr]))
  lmax = np.max(np.array([U.lp for U in U0_arr]))

  lmin = (1. - 0.01*np.sign(lmin)) * lmin
  lmax = (1. + 0.01*np.sign(lmin)) * lmax

  x0L = x - lmax * (teval - tstart)
  x0R = x - lmin * (teval - tstart)

  U0L = computeU0(U_ref, x0L)   # left  state at tstart
  U0R = computeU0(U_ref, x0R)   # right state at tstart

  x_tar  = -1   # target x
  xL_tar = -1
  xR_tar = -1

  # Running bisection
  c = 0
  while abs(x_tar - x) > x_prec:
    x0 = (x0L+x0R)/2.

    U0  = computeU0(U_ref, x0)
    U0L = computeU0(U_ref, x0L)
    U0R = computeU0(U_ref, x0R)

    x_tar  = x0  + U0.lp  * (teval - tstart)
    xL_tar = x0L + U0L.lp * (teval - tstart)
    xR_tar = x0R + U0R.lp * (teval - tstart)

    # checking that x is still between xL_tar and xR_tar:
    if (x-xL_tar)*(x-xR_tar) > 0.:
      print("x left the target interval! x = %e %i") %(x,c)
      exit(1)

    if x_tar < x:
      U0L = U0
      x0L = x0
    else:
      U0R = U0
      x0R = x0

    c += 1
    if c > 100:
      print("findU() timeout. Prec = %e") %(x_tar - x)
      break

  return(U0)



# ----------------------------------------------------------------------------------------
def computeL1(x, xth, dl):
  """
  Computes the L1 error on the specified quantity and returns it.
  Arguments:
  - x   : numerical solution
  - xth : exact solution 
  Returns:
  - L1  : The L1 error of the numerical solution
  """
  return np.sum(np.abs(x-xth) * dl) / np.sum(dl)
    # The dl comes from the fact that the L1 error is an integral and not a sum.
    # It reduces to a sum in the case of regular grid, but not for moving mesh.
    # (Fore more info: issue #21 on github)



# ----------------------------------------------------------------------------------------
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




# # ----------------------------------------------------------------------------------------
# # Script starts here
# # ----------------------------------------------------------------------------------------
# # args = parseArguments()

# # os.chdir(args.datadir)

# U_ref = State(rho_ref,v_ref,p_ref)
# U_ref.prim2aux()
# U_ref.Jm = Jm(U_ref)
# U_ref.Jp = Jp(U_ref)

# # Evaluating initial conditions for reference:
# # These values are used to constrain the bisection used to derive future fluid states.
# x0, U0 = fillU0(U_ref)

# # Solving for later conditions:
# # if not args.filename:
# x = setup_x()

# # else:
# #   file, status = openfile(args.filename)
# #   if status == 1:
# #     exit(1)

# #   time = file.get('/').attrs["time"]

# #   x   = file["/cells"]["x"]
# #   rho = file["/cells"]["rho"]
# #   v   = file["/cells"]["v"]
# #   p   = file["/cells"]["p"]
# #   D   = file["/cells"]["D"]
# #   E   = file["/cells"]["E"]
# #   dl  = file["/cells"]["dl"]

# # Computing exact solution:
# print("Computing exact solution:")
# U_arr = []
# for i in range(len(x)):
#   if i%10 == 0:
#     print("%i/%i") %(i, len(x))
#   U = findU(U_ref, x[i], time, U0)
#   U_arr.append(U)

# rhoth = np.array([U.rho for U in U_arr])
# vth   = np.array([U.v for U in U_arr])
# pth   = np.array([U.p for U in U_arr])


# # Computing the L1 error: (only possible when comparing with files)
# # L1 = computeL1(v,vth,dl)
# # print("t=%f, L1=%e") %(time, L1)


# pcolor = "tab:green"
# rcolor = "tab:orange"
# vcolor = "tab:blue"


# if args.plot:

#   if args.initial_state:
#     file_ini, status = openfile("phys0000000000.h5")
#     if status == 1:
#       exit(1)

#     time_ini = file_ini.get('/').attrs["time"]

#     x_ini   = file_ini["/cells"]["x"]
#     rho_ini = file_ini["/cells"]["rho"]
#     v_ini   = file_ini["/cells"]["v"]
#     p_ini   = file_ini["/cells"]["p"]
#     D_ini   = file_ini["/cells"]["D"]
#     E_ini   = file_ini["/cells"]["E"]
#     dl_ini  = file_ini["/cells"]["dl"]

#     rhoth_ini = np.array([U.rho for U in U0])
#     vth_ini   = np.array([U.v for U in U0])
#     pth_ini   = np.array([U.p for U in U0])

#     plt.scatter(x_ini,rho_ini/np.min(rhoth), marker='x', c=rcolor)
#     plt.plot(x0,rhoth_ini/np.min(rhoth), c='k')
#     plt.scatter(x_ini,v_ini, marker='+', c=vcolor)
#     plt.plot(x0,vth_ini, c='k')
#     plt.scatter(x_ini,p_ini/np.min(pth), marker='1', c=pcolor)
#     plt.plot(x0,pth_ini/np.min(pth), c='k')

#   # Plot:
#   plt.plot(x,rho/np.min(rhoth), label='rho, numerical', marker='x', c=rcolor)
#   plt.plot(x,rhoth/np.min(rhoth), c='k')
#   plt.plot(x,v, label='v, numerical', marker='+', c=vcolor)
#   plt.plot(x,vth, c='k')
#   plt.plot(x,p/np.min(pth), label='p, numerical', marker='1', c=pcolor)
#   plt.plot(x,pth/np.min(pth), c='k')

#   if args.cell_size:
#     plt.plot(x,dl/np.max(dl), label='dl/dlmax')

#   # plt.plot(x,(rho-rhoth)/np.min(rhoth), label='rho, residual', marker='+')
#   # plt.plot(x,(v-vth)/np.min(vth), label='v, residual', marker='+')
#   # plt.plot(x,(p-pth)/np.min(pth), label='p, residual', marker='+')

#   # plt.plot(x,np.zeros(len(x)),'k--')

#   plt.xlabel("x")
#   plt.title("Isentropic wave at t=0.8s")

#   plt.legend()
#   plt.show()



















