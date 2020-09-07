# -*- coding: utf-8 -*-
# @Author: eliotayache
# @Date:   2020-05-14 16:24:48
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-09-07 17:38:01


import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib.colors import LogNorm
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


def readData(key, it):
  data = pd.read_csv('../../results/temp/%s%d.out'  %(key,it), sep=" ")
  return(data)


def quadMesh(data, key, 
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

  plt.figure()
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

    if log==True:
      plt.pcolor(xj, yj, zj, 
        norm=LogNorm(vmin=vmin, vmax=vmax), 
        edgecolors=edges, 
        facecolor=color)
    else:
      plt.pcolor(xj, yj, zj, 
        vmin=vmin, vmax=vmax, 
        edgecolors=edges, 
        facecolor=color)

  plt.colorbar()

  if (quiver):
    if (geometry=='polar'):
      xx = x * np.cos(y)
      yy = x * np.sin(y)
      u = vx * np.cos(y)
      v = vx * np.sin(y)
      q = plt.quiver(xx,yy,u,v, alpha=0.2)

  plt.show()
