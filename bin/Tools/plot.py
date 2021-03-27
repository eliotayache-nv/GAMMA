# -*- coding: utf-8 -*-
# @Author: eliotayache
# @Date:   2020-05-14 16:24:48
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-03-27 22:11:58


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

from gamma_io import *
from isentropic import *
from BM import *
from synchrotron import *

# plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.rc('legend', fontsize=12) 
plt.rcParams['savefig.dpi'] = 300

def plot(key, it, contour=False):
  f = plt.figure()
  data = pd.read_csv('../../results/temp/%s%d.out' %(key,it), sep=" ", header=None)
  if contour:
    sns.set_style("white")
    sns.kdeplot(data)
  else:
    sns.heatmap(data.T); f.show()

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


def getArray(data, key):
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
  key2=None,
  log2=False,
  v1min=None,
  v2min=None,
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
  axis=None, 
  thetaobs=0.,
  nuobs=1.e17,
  expand=False):

  if key2 and geometry!="polar":
    print("Use polar geometry with key2")
    return

  if z_override is not None:
    z = z_override
  elif key == "PSyn":
    PSyn = addSynchPower(data, phi, thetaobs, nuobs)
    data["PSyn"] = PSyn
    z = data.pivot(index='j', columns='i', values="PSyn").to_numpy()
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

  # z = np.ma.masked_array(z, np.isnan(z))
  # x = np.ma.masked_array(x, np.isnan(x))
  # y = np.ma.masked_array(y, np.isnan(y))
  # dx = np.ma.masked_array(dx, np.isnan(dx))
  # dy = np.ma.masked_array(dy, np.isnan(dy))

  if key2:
    z2 = data.pivot(index='j', columns='i', values=key2).to_numpy()
    z2 = np.append(z2,np.expand_dims(z2[-1,:], axis=0), axis=0)
    z2 = np.ma.masked_array(z2, np.isnan(z2))

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
  if key2:
    vmin2 = np.nanmin(z2)
    vmax2 = np.nanmax(z2[4:,:])
  if v1min:
    vmin = v1min
  if v2min:
    vmin2 = v2min

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
    if key2:
      ax.set_thetamin(-ymax*180./np.pi)
    else:
      ax.set_thetamin(ymin*180./np.pi)

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

    # if invert==True:
    #   xj

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

    if key2:  
      if log2==True:
        im2 = ax.pcolor(-yj, xj, zj2, 
          norm=LogNorm(vmin=vmin2, vmax=vmax2), 
          edgecolors=edges,
          cmap=cmap,
          facecolor=color)
      else:
        im2 = ax.pcolor(-yj, xj, zj2, 
          vmin=vmin2, vmax=vmax2, 
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
  if geometry == "polar":
    # ax.set_aspect('equal')
    ax.set_rorigin(0)
    ax.set_rmin(xmin)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_rticks([xmin, xmax])

  if colorbar:
    cb = f.colorbar(im, ax=ax, orientation='vertical', shrink=.6, pad=0.1)
    cb.set_label(key, fontsize=14)
    if key2:
      cb2 = f.colorbar(im2, orientation='vertical')
      cb2.set_label(key2, fontsize=14)

  if tlayout:
    f.tight_layout()

  return xmin, xmax, im


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


# ----------------------------------------------------------------------------------------
# Specific functions
def isenwave(data):
  time = data["t"][0]
  f, axes = plotMulti(data, ["rho","p","vx"], 
    tracer=False, 
    line=False, 
    labels={"rho":"$\\rho$", "vx":"v"})
  plotIsen1D(time, "rho", ax=axes[0], color="r", label="exact")
  plotIsen1D(time, "p", ax=axes[1], color="r")
  plotIsen1D(time, "vx", ax=axes[2], color="r")
  plt.xlabel("x")
  axes[0].legend()
  plt.tight_layout()
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
def AnalyseBoxFit(data, jtrack=0, full=False):

  E0 = 1.e53
  n0 = 1.e0
  t = data["t"][0]
  BW = BM(E0, n0, t)
  RShock = BW.RShock/lNorm
  x_norm = None

  f, axes = plotMulti(data, ["rho","p","lfac",'gmax'], 
    jtrack=jtrack,
    tracer=False, 
    line=False, 
    labels={"rho":"$\\rho/\\rho_0$", 
            "p":"$p/p_0$",
            "lfac":"$\\gamma$",
            "gmax":"$\\gamma_\\{max}"}, 
    x_norm=x_norm)

  plotBM1D(data, "rho", jtrack=jtrack, x_norm=x_norm, ax=axes[0], color="r", label="BM", zorder=10)
  plotBM1D(data, "p", jtrack=jtrack, x_norm=x_norm, ax=axes[1], color="r", zorder=10)
  plotBM1D(data, "lfac", jtrack=jtrack, x_norm=x_norm, ax=axes[2], color="r", zorder=10)
  plotBM1D(data, "gmax", jtrack=jtrack, x_norm=x_norm, ax=axes[3], color="r", zorder=10)

  plt.xlim(0.90*RShock, 1.05*RShock)
  axes[0].set_yscale("log")
  axes[1].set_yscale("log")
  axes[2].set_yscale("log")
  axes[3].set_yscale("log")
  plt.xlabel("$r [cm]$")
  axes[0].legend()
  plt.tight_layout()
  
  if full==True:
    return(BoxFitImages(data))


def BoxFitImages(data):

  import matplotlib.ticker as ticker
  frmtr = ticker.FormatStrFormatter('%4.1e')

  thetamax = 0.5
  rmin, rmax, im = quadMesh(data, "rho", log=True, geometry='polar')
  rp = (rmax - rmin*np.cos(thetamax))
  h = rmax*np.sin(thetamax)
  fig, [ax1, ax2] = plt.subplots(2, 1, 
                                 subplot_kw=dict(polar=True), 
                                 figsize=(rp / rmin*3 * 1.3, 2*h / rmin*3))
  rmin, rmax, im1 = quadMesh(data, "rho", log=True, fig=fig, axis=ax1, colorbar=False)
  cb = fig.colorbar(im1, ax=ax1, orientation='vertical', shrink=.7, pad=0.1)
  cb.set_label("density", fontsize=14)
  cax1 = cb.ax
  ax1.set_thetalim(0, thetamax)
  ax1.set_rorigin(0)
  ax1.set_rmin(rmin)
  ax1.set_rmax(rmax)
  ax1.tick_params(labelleft=False, labelright=True, labeltop=False, labelbottom=True)
  ax1.yaxis.set_major_formatter(frmtr)

  rmin, rmax, im2 = quadMesh(data, "p", log=True, fig=fig, axis=ax2, cmap="cividis", colorbar=False)
  cb = fig.colorbar(im2, ax=ax2, orientation='vertical', shrink=.7, pad=0.1, cmap="cividis")
  cb.set_label("pressure", fontsize=14)
  cax2 = cb.ax
  ax2.set_theta_direction(-1)
  ax2.set_thetalim(0, thetamax)
  ax2.set_rorigin(0)
  ax2.set_rmin(rmin)
  ax2.set_rmax(rmax)
  ax2.tick_params(labelleft=False, labelright=True, labeltop=False, labelbottom=True)
  ax2.yaxis.set_major_formatter(frmtr)

  ax2.set_xticks(ax1.get_xticks()[1:])
  plt.subplots_adjust(hspace=0, bottom=0.)

  s = 2./rp
  d = 0.9*(2-h*s)/2.

  pos1 = ax1.get_position()
  pos2 = ax2.get_position()
  cpos1 = cax1.get_position()
  cpos2 = cax2.get_position()
  ratio1 = (pos1.height)/2.
  ratio2 = (pos2.height)/2.
  ax1.set_position([pos1.x0, pos1.y0-d*ratio1, pos1.width, pos1.height])
  ax2.set_position([pos2.x0, pos2.y0+d*ratio2, pos2.width, pos2.height])
  cax1.set_position([cpos1.x0, cpos1.y0-d*ratio1, cpos1.width, cpos1.height])
  cax2.set_position([cpos2.x0, cpos2.y0+d*ratio2, cpos2.width, cpos2.height])
    # /4 because ratio of 0.5/2

  thetamax = 0.5
  rmin, rmax, im = quadMesh(data, "gmax", log=True, geometry='polar')
  rp = (rmax - rmin*np.cos(thetamax))
  h = rmax*np.sin(thetamax)
  fig, [ax1, ax2] = plt.subplots(2, 1, 
                                 subplot_kw=dict(polar=True), 
                                 figsize=(rp / rmin*3 * 1.3, 2*h / rmin*3))
  rmin, rmax, im1 = quadMesh(data, "gmax", log=True, fig=fig, axis=ax1, colorbar=False)
  cb = fig.colorbar(im1, ax=ax1, orientation='vertical', shrink=.7, pad=0.1)
  cb.set_label("$\\gamma_\\mathrm{max}$", fontsize=14)
  cax1 = cb.ax
  ax1.set_thetalim(0, thetamax)
  ax1.set_rorigin(0)
  ax1.set_rmin(rmin)
  ax1.set_rmax(rmax)
  ax1.tick_params(labelleft=False, labelright=True, labeltop=False, labelbottom=True)
  ax1.yaxis.set_major_formatter(frmtr)

  rmin, rmax, im2 = quadMesh(data, "psyn", v1min=2, fig=fig, axis=ax2, cmap="cividis", colorbar=False)
  cb = fig.colorbar(im2, ax=ax2, orientation='vertical', shrink=.7, pad=0.1, cmap="cividis")
  cb.set_label("spectral index p", fontsize=14)
  cax2 = cb.ax
  ax2.set_theta_direction(-1)
  ax2.set_thetalim(0, thetamax)
  ax2.set_rorigin(0)
  ax2.set_rmin(rmin)
  ax2.set_rmax(rmax)
  ax2.tick_params(labelleft=False, labelright=True, labeltop=False, labelbottom=True)
  ax2.yaxis.set_major_formatter(frmtr)

  ax2.set_xticks(ax1.get_xticks()[1:])
  plt.subplots_adjust(hspace=0, bottom=0.)

  s = 2./rp
  d = 0.9*(2-h*s)/2.

  pos1 = ax1.get_position()
  pos2 = ax2.get_position()
  cpos1 = cax1.get_position()
  cpos2 = cax2.get_position()
  ratio1 = (pos1.height)/2.
  ratio2 = (pos2.height)/2.
  ax1.set_position([pos1.x0, pos1.y0-d*ratio1, pos1.width, pos1.height])
  ax2.set_position([pos2.x0, pos2.y0+d*ratio2, pos2.width, pos2.height])
  cax1.set_position([cpos1.x0, cpos1.y0-d*ratio1, cpos1.width, cpos1.height])
  cax2.set_position([cpos2.x0, cpos2.y0+d*ratio2, cpos2.width, cpos2.height])
    # /4 because ratio of 0.5/2
  




























