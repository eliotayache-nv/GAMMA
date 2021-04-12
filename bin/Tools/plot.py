# -*- coding: utf-8 -*-
# @Author: eliotayache
# @Date:   2020-05-14 16:24:48
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-04-12 14:24:40

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import glob
import os
import string
from astropy.io import ascii

from gamma_io import *
from isentropic import *
from BM import *
from synchrotron import *

# plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.rc('legend', fontsize=12) 
plt.rcParams['savefig.dpi'] = 200

# def plot(key, it, contour=False):
#   f = plt.figure()
#   data = pd.read_csv('../../results/temp/%s%d.out' %(key,it), sep=" ", header=None)
#   if contour:
#     sns.set_style("white")
#     sns.kdeplot(data)
#   else:
#     sns.heatmap(data.T); f.show()

def tricontour(name, it, key,):
  data = pd.read_csv('../../results/temp/%s%d.out'  %(name,it), sep=" ")

  x = data['x']
  y = data['y']
  z = np.log10(data[key])

  fig = plt.figure()
  ax = plt.gca()

  ax.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
  cntr = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

  fig.colorbar(cntr, ax=ax)
  ax.plot(x, y, 'ko', ms=0.3)
  ax.set(xlim=(-2, 2), ylim=(-2, 2))

  return x,y,z


def getArray(data, key, oneDimensional=False):
  if oneDimensional:
    return(data[key].to_numpy())
  else:
    return(data.pivot(index='j', columns='i', values=key).to_numpy())



def plotMulti(data, keys, jtrack=None, log=[], labels={}, **kwargs):

  Nk = len(keys)
  f, axes = plt.subplots(Nk, 1, sharex=True, figsize=(6,2*Nk))

  for key,k,ax in zip(keys,range(Nk),axes):

    logkey = False
    label = None
    if key in log:
      logkey = True
    if key in labels:
      label = labels[key]

    plot1D(data, key, ax, jtrack=jtrack, log=logkey, label=label, **kwargs)

  plt.tight_layout()

  return(f, axes)



def plot1D(data, key, 
  ax = None,
  mov="x",
  log=False, 
  v1min=None,
  tracer=True,
  line=True,
  r2=False,
  x_norm=None,
  jtrack=None,
  label=None, **kwargs):

  if key=="lfac":
    var = "vx"
  else:
    var = key

  if jtrack is not(None):
    z = pivot(data, var)[jtrack,:]
    x = pivot(data, "x")[jtrack,:]
    tracvals = pivot(data, "trac")[jtrack,:]
  else:
    z = data[var].to_numpy()
    x  = np.copy(data["x"].to_numpy())
    tracvals = data["trac"].to_numpy()

  if x_norm is not(None):
    x /= x_norm

  if key=="lfac":
    z = 1./np.sqrt(1 - z**2)

  if r2:
    z*=x**2

  xmin = np.min(x)
  xmax = np.max(x)

  vmin = np.min(z)
  vmax = np.max(z)

  if v1min:
    vmin = v1min

  if ax == None:
    f = plt.figure()
    ax = plt.gca()

  if label is not(None):
    ax.set_ylabel(label)
  else:
    ax.set_ylabel(key)

  if log==True:
    ax.set_yscale('log')

  if line:
    ax.plot(x,z,'k',zorder=1)    
  ax.scatter(x,z, c='None', edgecolors='k', lw=2, zorder=2, label="numerical")
  if tracer:
    ax.scatter(x,z, c=tracvals, edgecolors='None', zorder=3, cmap='cividis')



def quadMesh(data, key,
  z_override = None,
  mov="x",
  log=False, 
  v1min=None,
  geometry="cartesian", 
  quiver=False, 
  color=None, 
  edges='None',
  invert=False,
  r2=False,
  cmap='magma',
  tlayout=False,
  colorbar=True,
  slick=False,
  phi=0.,
  fig=None,
  label=None,
  axis=None, 
  thetaobs=0.,
  nuobs=1.e17,
  shrink=0.6,
  expand=False):

  # if key2 and geometry!="polar":
  #   print("Use polar geometry with key2")
  #   return

  if z_override is not None:
    z = z_override
  elif key == "PSyn":
    PSyn = addSynchPower(data, phi, thetaobs, nuobs)
    data["PSyn"] = PSyn
    z = data.pivot(index='j', columns='i', values="PSyn").to_numpy()
    zmin = np.min(z[np.logical_and(~np.isnan(z), z>0)])
    z[np.logical_or(np.isnan(z), z<zmin)] = zmin
  elif key == "lfac":
    vx = data.pivot(index='j', columns='i', values="vx").to_numpy()
    vy = data.pivot(index='j', columns='i', values="vy").to_numpy()
    z = 1./np.sqrt(1 - (vx**2+vy**2))    
  else:
    z = data.pivot(index='j', columns='i', values=key).to_numpy()
  x  = data.pivot(index='j', columns='i', values='x').to_numpy()
  dx = data.pivot(index='j', columns='i', values='dx').to_numpy()
  y  = data.pivot(index='j', columns='i', values='y').to_numpy()
  dy = data.pivot(index='j', columns='i', values='dy').to_numpy()

  # duplicating last row for plotting
  z = np.append(z,np.expand_dims(z[-1,:], axis=0), axis=0)
  x = np.append(x,np.expand_dims(x[-1,:], axis=0), axis=0)
  dx = np.append(dx,np.expand_dims(dx[-1,:], axis=0), axis=0)
  y = np.append(y,np.expand_dims(y[-1,:], axis=0), axis=0)
  dy = np.append(dy,np.expand_dims(dy[-1,:], axis=0), axis=0)

  # duplicating first column for plotting
  z = np.append(z, np.expand_dims(z[:,-1], axis=1), axis=1)
  x = np.append(x, np.expand_dims(x[:,-1], axis=1), axis=1)
  dx = np.append(dx, np.expand_dims(dx[:,-1], axis=1), axis=1)
  y = np.append(y, np.expand_dims(y[:,-1], axis=1), axis=1)
  dy = np.append(dy, np.expand_dims(dy[:,-1], axis=1), axis=1)

  nact = np.array([np.count_nonzero(~np.isnan(xj)) for xj in x])

  if (quiver):
    vx = data.pivot(index='j', columns='i', values='vx').to_numpy()
    vx = np.append(vx,np.expand_dims(vx[-1,:], axis=0), axis=0)
    vx = np.ma.masked_array(vx, np.isnan(vx))

  if r2:
    z*=x**2

  xmin = np.nanmin(x)
  xmax = np.nanmax(x)
  ymin = np.nanmin(y)
  ymax = np.nanmax(y)

  vmax = np.nanmax(z[4:,:])
  vmin = np.nanmin(z)
  if log == True:
    vmin = np.nanmin(z[z>0])
  if v1min:
    vmin = v1min

  if geometry == "polar":
    projection="polar"
  else:
    projection=None

  if axis is None:
    f = plt.figure()
    ax = plt.axes(projection=projection)
  else:
    f = fig
    ax = axis

  if geometry == "polar" or axis is not None:
    ax.set_thetamax(ymax*180./np.pi)
    ax.set_thetamin(ymin*180./np.pi)
    if invert==True:
      ax.set_thetamin(-ymax*180./np.pi)

  if slick:
    ax.axis("off")

  for j in range(z.shape[0]-1):
    xj = x - dx/2.
    yj = y - dy/2.
    dyj = dy
    xj[j,nact[j]-1] += dx[j,nact[j]-1]

    if mov == 'y':
      tmp = np.copy(xj)
      xj = yj
      yj = np.copy(tmp)
      dyj = dx

    xj[j+1,:] = xj[j,:]   
    yj[j+1,:] = yj[j,:]+dyj[j,:]   

    xj = xj[j:j+2,:]
    yj = yj[j:j+2,:]
    zj = z[j:j+2,:]

    if invert==True:
      yj *= -1

    if log==True:
      im = ax.pcolor(yj, xj, zj, 
        norm=LogNorm(vmin=vmin, vmax=vmax), 
        edgecolors=edges,
        cmap=cmap,
        facecolor=color)
    else:
      im = ax.pcolor(yj, xj, zj, 
        vmin=vmin, vmax=vmax, 
        edgecolors=edges, 
        cmap=cmap,
        facecolor=color)

  # if (quiver):
    # if (geometry=='polar'):
    #   xx = x * np.cos(y)
    #   yy = x * np.sin(y)
    #   u = vx * np.cos(y)
    #   v = vx * np.sin(y)
    #   q = plt.quiver(yy[:,::5],xx[:,::5],v[:,::5],u[:,::5], headwidth=10, headlength=10)

  if geometry != "polar":
    ax.set_aspect('equal')
  if geometry == "polar" or axis is not None:
    ax.set_rorigin(0)
    ax.set_rmin(xmin)
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_rticks([xmin, xmax])

  if colorbar:
    cb = f.colorbar(im, ax=ax, orientation='vertical', shrink=shrink, pad=0.1)
    if label is None:
      label = key  
    cb.set_label(label, fontsize=14)

  if tlayout:
    f.tight_layout()

  thetamax = ymax*180./np.pi
  return xmin, xmax, thetamax, im


def loopFigs(dir, key, oneDimensional = False, **kwargs):
  if not os.path.exists("../../results/%s/figs" %dir):
    os.mkdir("../../results/%s/figs" %dir)
  for filename in glob.glob("../../results/%s/*.out" %dir):
    print(filename.split("results/")[1])
    data = readData(filename.split("results/")[1])
    if oneDimensional:
      plot1D(data, key, **kwargs)      
    else:
      quadMesh(data, key, **kwargs)
    plt.savefig("../../results/%s/figs/%s.png" %(dir,filename.rstrip(".out").split(dir)[1]))
    plt.close()


def plotConvergence(res,err):
  res = np.array(res)
  err = np.array(err)
  plt.plot(res, err, marker='o', mec='k', ls='--', c='k')
  slopes = np.gradient(np.log10(err), np.log10(res))
  print(slopes)
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('resolution\n (initial no of cells)')
  plt.ylabel('L1 error')



# ----------------------------------------------------------------------------------------
# Specific functions
def isenwave(data):
  time = data["t"][0]
  f, axes = plotMulti(data, ["rho","p","vx"], 
    tracer=False, 
    line=False, 
    labels={"rho":"density",'p':'pressure', "vx":"velocity"})
  plotIsen1D(time, "rho", ax=axes[0], color="r", label="exact")
  plotIsen1D(time, "p", ax=axes[1], color="r")
  plotIsen1D(time, "vx", ax=axes[2], color="r")
  plt.xlabel("x")
  axes[0].legend()
  plt.tight_layout()

  computeOrderOfPrec(data)

  return(f, axes)


# ----------------------------------------------------------------------------------------
# Specific functions
def BMwave(data):

  E0 = 1.e53
  n0 = 1.e0
  t = data["t"][0]
  BW = BM(E0, n0, t)
  RShock = BW.RShock/lNorm

  f, axes = plotMulti(data, ["rho","p","lfac"], 
    tracer=False, 
    line=False, 
    labels={"rho":"$\\rho/\\rho_0$", "p":"$p/p_0$","lfac":"$\\gamma$"}, 
    x_norm=RShock)

  plotBM1D(data, "rho", x_norm=RShock, ax=axes[0], color="r", label="exact", zorder=10)
  plotBM1D(data, "p", x_norm=RShock, ax=axes[1], color="r", zorder=10)
  plotBM1D(data, "lfac", x_norm=RShock, ax=axes[2], color="r", zorder=10)

  plt.xlim(0.997, 1.001)
  axes[0].set_yscale("log")
  axes[1].set_yscale("log")
  axes[2].set_yscale("log")
  plt.xlabel("$r/r_\\mathrm{shock}$")
  axes[0].legend()
  plt.tight_layout()
  return(f, axes)


# ----------------------------------------------------------------------------------------
# Specific functions
def AnalyseBoxFit(data, jtrack=0, full=False, thetamax=0.5):

  E0 = 1.e53
  n0 = 1.e0
  t = data["t"][0]
  BW = BM(E0, n0, t)
  RShock = BW.RShock/lNorm
  x_norm = None

  f, axes = plotMulti(data, ["rho","p","lfac",'gmax','gmin'], 
    jtrack=jtrack,
    tracer=False, 
    line=False, 
    # labels={"rho":"$\\rho/\\rho_0$", 
    #         "p":"$p/p_0$",
    #         "lfac":"$\\gamma$",
    #         "gmax":"$\\gamma_\\mathrm{max}$",
    #         "gmin":"$\\gamma_\\mathrm{min}$"}, 
    labels={"rho":"density\n(normalised)", 
            "p":"pressure\n(normalised)",
            "lfac":"Lorentz factor",
            "gmax":"$\\gamma_\\mathrm{max}$",
            "gmin":"$\\gamma_\\mathrm{min}$"}, 
    x_norm=x_norm)

  plotBM1D(data, "rho", jtrack=jtrack, x_norm=x_norm, ax=axes[0], color="r", label="BM", zorder=10)
  plotBM1D(data, "p", jtrack=jtrack, x_norm=x_norm, ax=axes[1], color="r", zorder=10)
  plotBM1D(data, "lfac", jtrack=jtrack, x_norm=x_norm, ax=axes[2], color="r", zorder=10)
  plotBM1D(data, "gmax", jtrack=jtrack, x_norm=x_norm, ax=axes[3], color="r", zorder=10)
  plotBM1D(data, "gmin", jtrack=jtrack, x_norm=x_norm, ax=axes[4], color="r", zorder=10)

  plt.xlim(0.90*RShock, 1.05*RShock)
  axes[0].set_yscale("log")
  axes[1].set_yscale("log")
  axes[2].set_yscale("log")
  axes[3].set_yscale("log")
  axes[4].set_yscale("log")
  plt.xlabel("$r$ (light seconds)")
  axes[0].legend()
  plt.tight_layout()
  
  if full==True:
    BoxFitImages(data, thetamax=thetamax)






def BoxFitImages(data, save=False, thetamax=0.5, rmin0=None, rmax0=None, **kwargs):

  thetamax *= 180./np.pi

  f = plt.figure()
  ax2 = plt.axes(projection='polar', frameon=False)
  ax = plt.axes(projection='polar')
  rmin, rmax, thetamaxim, im1 = quadMesh(data, "rho", log=True,
                                       fig=f, axis=ax, colorbar=False, **kwargs)
  rmin, rmax, thetamaxim, im2 = quadMesh(data, "p", log=True,
                                       fig=f, axis=ax, 
                                       invert=True, cmap='cividis', colorbar=False, **kwargs)
  ax.axvline(0, color='k', lw=0.7)

  if rmin0 is not None:
    rmin = rmin0
  if rmax0 is not None:
    rmax = rmax0

  ax.set_rmin(rmin)
  rmid = (rmax+rmin)/2. 
  ax.set_rticks([rmin, rmid, rmax])
  ax.set_thetamin(-thetamax)
  ax.set_thetamax(thetamax)

  ax2.set_rorigin(0)
  ax2.set_rmin(rmin)
  ax2.set_rmax(rmax)
  ax2.set_thetamin(-thetamax)
  ax2.set_thetamax(thetamax)

  cb = f.colorbar(im1, ax=ax, orientation='vertical', pad=0.1)
  cb.set_label("density", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.52, pos.width, 0.23])

  cb = f.colorbar(im2, ax=ax2, orientation='vertical', pad=0.1)
  cb.set_label("pressure", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.25, pos.width, 0.23])

  time = data['t'][0]
  plt.title('-\ntime = %.2e s\n radius in cm' %time, fontsize=12)

  if save==True:
    f.savefig('boxfit.png', bbox_inches='tight')


  f = plt.figure()
  ax2 = plt.axes(projection='polar', frameon=False)
  ax = plt.axes(projection='polar')
  rmin, rmax, thetamaxim, im1 = quadMesh(data, "gmax", log=True,
                                       fig=f, axis=ax, colorbar=False)
  rmin, rmax, thetamaxim, im2 = quadMesh(data, "gmin", log=True, 
                                       fig=f, axis=ax, 
                                       invert=True, cmap='cividis', colorbar=False)
  ax.axvline(0, color='k', lw=0.7)

  ax.set_rmin(rmin)
  rmid = (rmax+rmin)/2. 
  ax.set_rticks([rmin, rmid, rmax])
  ax.set_thetamin(-thetamax)
  ax.set_thetamax(thetamax)

  ax2.set_rorigin(0)
  ax2.set_rmin(rmin)
  ax2.set_rmax(rmax)
  ax2.set_thetamin(-thetamax)
  ax2.set_thetamax(thetamax)

  cb = f.colorbar(im1, ax=ax, orientation='vertical', pad=0.1)
  cb.set_label("$\\gamma_\\mathrm{max}$", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.52, pos.width, 0.23])

  cb = f.colorbar(im2, ax=ax2, orientation='vertical', pad=0.1)
  cb.set_label("$\\gamma_\\mathrm{min}$", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.25, pos.width, 0.23])

  plt.title('-\ntime = %.2e s\n radius in cm' %time, loc='center', fontsize=12)

  if save==True:
    f.savefig('test.png', bbox_inches='tight')




  f = plt.figure()
  ax2 = plt.axes(projection='polar', frameon=False)
  ax = plt.axes(projection='polar')
  rmin, rmax, thetamaxim, im1 = quadMesh(data, "PSyn", log=True,
                                       fig=f, axis=ax, colorbar=False)
  rmin, rmax, thetamaxim, im2 = quadMesh(data, "psyn", v1min=2., 
                                       fig=f, axis=ax, 
                                       invert=True, cmap='cividis', colorbar=False)
  ax.axvline(0, color='k', lw=0.7)

  ax.set_rmin(rmin)
  rmid = (rmax+rmin)/2. 
  ax.set_rticks([rmin, rmid, rmax])
  ax.set_thetamin(-thetamax)
  ax.set_thetamax(thetamax)

  ax2.set_rorigin(0)
  ax2.set_rmin(rmin)
  ax2.set_rmax(rmax)
  ax2.set_thetamin(-thetamax)
  ax2.set_thetamax(thetamax)

  cb = f.colorbar(im1, ax=ax, orientation='vertical', pad=0.1)
  cb.set_label("X-ray Emissivity\n (code units)", fontsize=10)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.52, pos.width, 0.23])

  cb = f.colorbar(im2, ax=ax2, orientation='vertical', pad=0.1)
  cb.set_label("spectral index", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.25, pos.width, 0.23])

  plt.title('-\ntime = %.2e s\n radius in cm' %time, loc='center', fontsize=12)

  if save==True:
    f.savefig('test.png', bbox_inches='tight')








# ----------------------------------------------------------------------------------------
# Specific functions
def AnalyseRT(data, save=False):

  f = plt.figure()
  ax2 = plt.axes(projection='polar', frameon=False)
  ax = plt.axes(projection='polar')
  rmin, rmax, thetamax, im1 = quadMesh(data, "rho", 
                                       log=True, fig=f, axis=ax, colorbar=False)
  rmin, rmax, thetamax, im2 = quadMesh(data, "trac", 
                                       fig=f, axis=ax, 
                                       invert=True, cmap='cividis', colorbar=False)
  ax.axvline(0, color='k', lw=0.7)

  rmin = (rmax+rmin)/2. 
  ax.set_rmin(rmin)
  ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
  rmid = (rmax+rmin)/2. 
  ax.set_rticks([rmin, rmid, rmax])

  ax2.set_rorigin(0)
  ax2.set_rmin(rmin)
  ax2.set_rmax(rmax)
  ax2.set_thetalim(-thetamax, thetamax)

  cb = f.colorbar(im1, ax=ax, orientation='vertical', pad=0.1)
  cb.set_label("density", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.52, pos.width, 0.23])

  cb = f.colorbar(im2, ax=ax2, orientation='vertical', pad=0.1)
  cb.set_label("tracer", fontsize=12)
  cax = cb.ax
  pos = cax.get_position()
  cax.set_position([pos.x0, 0.25, pos.width, 0.23])

  ax.set_xlabel('radius (cm)')
  time = data['t'][0]
  plt.title('time = %.2e s' %time, x=0.2, y=0.8, fontsize=12)

  if save==True:
    f.savefig('test.png', bbox_inches='tight')






















