# -*- coding: utf-8 -*-
# @Author: eliotayache
# @Date:   2020-05-14 16:24:48
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-09-08 12:18:07


import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import glob
import os

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


def readData(key, it, sequence=False):
  if sequence:
    filename = '../../results/%s%010d.out'  %(key,it)
  else:
    filename = '../../results/%s%d.out'  %(key,it)
  data = pd.read_csv(filename, sep=" ")
  return(data)


def quadMesh(data, key, 
  key2=None,
  log2=False,
  geometry="cartesian", 
  quiver=False, 
  color=None, 
  log=False, 
  edges='None',
  r2=False):

  z = data.pivot(index='j', columns='i', values=key).to_numpy()
  x  = data.pivot(index='j', columns='i', values='x').to_numpy()
  dx = data.pivot(index='j', columns='i', values='dx').to_numpy()
  y  = data.pivot(index='j', columns='i', values='y').to_numpy()
  dy = data.pivot(index='j', columns='i', values='dy').to_numpy()

  if key2:
    z2 = data.pivot(index='j', columns='i', values=key2).to_numpy()

  if (quiver):
    vx = data.pivot(index='j', columns='i', values='vx').to_numpy()

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

  if (key2):
    f, (ax2, ax) = plt.subplots(1,2, sharey=True)
  else:
    f = plt.figure()
    ax = plt.gca()

  for j in range(z.shape[0]-1):
    xj = x - dx/2.
    yj = y - dy/2.
    xj[j+1,:] = xj[j,:]
    
    if (geometry=='polar'):
      xx = xj * np.cos(yj)
      yy = xj * np.sin(yj)
      xj = xx
      yj = yy

    mask = np.zeros(z.shape)+1
    mask[j,:] = 0
    zj = np.ma.masked_array(z, mask>0)
    zj2 = np.ma.masked_array(z2, mask>0)

    if log==True:
      im = ax.pcolor(yj, xj, zj, 
        norm=LogNorm(vmin=vmin, vmax=vmax), 
        edgecolors=edges,
        facecolor=color)
    else:
      im = ax.pcolor(yj, xj, zj, 
        vmin=vmin, vmax=vmax, 
        edgecolors=edges, 
        facecolor=color)

    if key2:  
      if log2==True:
        im2 = ax2.pcolor(-yj, xj, zj2, 
          norm=LogNorm(vmin=vmin2, vmax=vmax2), 
          edgecolors=edges,
          facecolor=color)
      else:
        im2 = ax2.pcolor(-yj, xj, zj2, 
          vmin=vmin2, vmax=vmax2, 
          edgecolors=edges, 
          facecolor=color)


  if (quiver):
    if (geometry=='polar'):
      xx = x * np.cos(y)
      yy = x * np.sin(y)
      u = vx * np.cos(y)
      v = vx * np.sin(y)
      q = plt.quiver(xx,yy,u,v, alpha=0.2)

  ax.set_aspect('equal')
  ax2.set_aspect('equal')
  f.colorbar(im, ax=ax)
  # f.colorbar(im2, ax=ax2)
  cb = plt.colorbar(im2,ax=[ax2],location='left')
  # f.tight_layout()


def loopFigs(dir, key, **kwargs):
  for filename in os.listdir("../../results/%s" %dir):
    print("../../results/%s/%s" %(dir,filename))
    data = readData("%s/%d" %(dir,filename))
    quadMesh(data)


