# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2020-12-02 15:00:33
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2020-12-02 17:47:08


import numpy as np
import pandas as pd


def readData(key, it=None, sequence=False):
  if sequence:
    filename = '../../results/%s/phys%010d.out'  %(key,it)
  elif it==None:
    filename = '../../results/%s'  %(key)
  else:
    filename = '../../results/%s%d.out'  %(key,it)
  data = pd.read_csv(filename, sep=" ")
  return(data)


def to_2d(folder, it, ny, ymin, ymax, geometry="cartesian"):

  data = readData(folder, it, True)
  data_2d = pd.DataFrame(np.tile(data.values,[ny,1]))
  data_2d.columns = data.columns
  nx = data.shape[0]
  
  j = np.arange(ny)
  jcol = np.repeat(j, nx)

  dy = (ymax-ymin)/ny
  dycol = np.repeat(dy, nx*ny)

  y = np.linspace(ymin+dy/2., ymax-dy/2., ny)
  ycol = np.repeat(y, nx)

  if geometry=="cartesian":
    dlycol = dycol
  if geometry=="polar":
    x = data_2d.x.to_numpy()
    dlycol = dycol*x

  data_2d.j = jcol
  data_2d.y = ycol
  data_2d.dy = dycol
  data_2d.dly = dlycol

  return(data_2d)


