# -*- coding: utf-8 -*-
# @Author: eliotayache
# @Date:   2020-05-14 16:24:48
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-11-23 21:30:15


import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import glob
import os
import string
import glob

# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.rcParams['savefig.dpi'] = 200


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


def readData(key, it=None, sequence=False):
  if sequence:
    filename = '../../results/%s/phys%010d.out'  %(key,it)
  elif it==None:
    filename = '../../results/%s'  %(key)
  else:
    filename = '../../results/%s%d.out'  %(key,it)
  data = pd.read_csv(filename, sep=" ")
  return(data)


def getArray(data, key):
  return(data.pivot(index='j', columns='i', values=key).to_numpy())


def plot1D(data, key, 
  mov="x",
  log=False, 
  v1min=None,
  color=None, 
  r2=False):

  if key=="lfac":
    var = "vx"
  else:
    var = key

  z = data[var].to_numpy()
  x  = data["x"].to_numpy()

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

  f = plt.figure()

  if log==True:
    plt.yscale('log')

   plt.plot(x,z,'ko-', mfc='r', mec='r', ms=3)



def quadMesh(data, key, 
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
  r2=False,
  cmap='magma',
  tlayout=True,
  colorbar=True):

  z = data.pivot(index='j', columns='i', values=key).to_numpy()
  x  = data.pivot(index='j', columns='i', values='x').to_numpy()
  dx = data.pivot(index='j', columns='i', values='dx').to_numpy()
  y  = data.pivot(index='j', columns='i', values='y').to_numpy()
  dy = data.pivot(index='j', columns='i', values='dy').to_numpy()
  z = np.ma.masked_array(z, np.isnan(z))
  x = np.ma.masked_array(x, np.isnan(x))
  y = np.ma.masked_array(y, np.isnan(y))
  dx = np.ma.masked_array(dx, np.isnan(dx))
  dy = np.ma.masked_array(dy, np.isnan(dy))

  if key2:
    z2 = data.pivot(index='j', columns='i', values=key2).to_numpy()
    z2 = np.ma.masked_array(z2, np.isnan(z2))

  if (quiver):
    vx = data.pivot(index='j', columns='i', values='vx').to_numpy()
    vx = np.ma.masked_array(vx, np.isnan(vx))

  if r2:
    z*=x**2

  xmin = np.min(x)
  xmax = np.max(x)
  ymin = np.min(y)
  ymax = np.max(y)

  vmin = np.min(z)
  vmax = np.max(z[4:,:])
  if key2:
    vmin2 = np.min(z2)
    vmax2 = np.max(z2[4:,:])
  if v1min:
    vmin = v1min
  if v2min:
    vmin2 = v2min

  if (key2):
    f, (ax2, ax) = plt.subplots(1,2, sharey=True)
  else:
    f = plt.figure()
    ax = plt.gca()

  for j in range(z.shape[0]-1):
    xj = x - dx/2.
    yj = y - dy/2.
    dyj = dy

    if mov == 'y':
      tmp = np.copy(xj)
      xj = yj
      yj = np.copy(tmp)
      dyj = dx

    xj[j+1,:] = xj[j,:]   #these xj[] might not exist
    yj[j+1,:] = yj[j,:]+dyj[j,:]   #these xj[] might not exist

    if (geometry=='polar'):
      xx = xj * np.cos(yj)
      yy = xj * np.sin(yj)
      xj = xx
      yj = yy


    mask = np.zeros(z.shape)+1
    mask[j,:] = 0
    mask[np.isnan(z)]=1
    zj = np.ma.masked_array(z, mask>0)
    if key2:
      zj2 = np.ma.masked_array(z2, mask>0)

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
        im2 = ax2.pcolor(-yj, xj, zj2, 
          norm=LogNorm(vmin=vmin2, vmax=vmax2), 
          edgecolors=edges,
          cmap=cmap,
          facecolor=color)
      else:
        im2 = ax2.pcolor(-yj, xj, zj2, 
          vmin=vmin2, vmax=vmax2, 
          edgecolors=edges, 
          cmap=cmap,
          facecolor=color)


  if (quiver):
    if (geometry=='polar'):
      xx = x * np.cos(y)
      yy = x * np.sin(y)
      u = vx * np.cos(y)
      v = vx * np.sin(y)
      q = plt.quiver(yy[:,::5],xx[:,::5],v[:,::5],u[:,::5], headwidth=10, headlength=10)

  ax.set_aspect('equal')
  if key2:
    ax2.set_aspect('equal')

  if colorbar:
    cb = f.colorbar(im, ax=ax, orientation='horizontal')
    cb.set_label(key, fontsize=18)
    if key2:
      cb2 = f.colorbar(im2, ax=ax2, orientation='horizontal')
      cb2.set_label(key2, fontsize=18)

  if tlayout:
    f.tight_layout()


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




