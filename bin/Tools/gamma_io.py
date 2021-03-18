# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2021-01-30 17:59:32
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-03-18 17:14:54

import pandas as pd
import numpy as np
import h5py
import os
import glob

mp_ = 1.6726219e-24
c_  = 2.99792458e10


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


def extractDigits(s, leading_zeros=True):
  if leading_zeros:
    return(''.join(filter(str.isdigit, s)))
  else:
    res = ''.join(filter(str.isdigit, s)).lstrip('0')
    if res != '':
      return res
    else:
      return '0'


def readData(key, it=None, sequence=False):
  if sequence:
    filename = '../../results/%s/phys%010d.out'  %(key,it)
  elif it==None:
    filename = '../../results/%s'  %(key)
  else:
    filename = '../../results/%s%d.out'  %(key,it)
  data = pd.read_csv(filename, sep=" ")
  return(data)


def applyToAll(dir, func, *args, **kwargs):
  filelist = sorted(glob.glob("../../results/%s/*.out" %dir))
  l=[]
  for filename,index in zip(filelist,range(len(filelist))):
    name = filename.split("/")[-1]
    print("%d, %s" %(index,name))
    it = int(extractDigits(name, leading_zeros=False))
    l.append(func(dir, it, index, *args, **kwargs))
  l = np.array(l)
  return(filelist, l) 


def extract(data, key):
  return(data[key].to_numpy())



def pivot(data, key):
  return(data.pivot(index="j", columns="i", values=key).to_numpy())



def pandas2double(data):

  T = data["t"][0] # time

  rho = extract(data, "rho")*mp_
  p = extract(data, "p")*mp_*c_**2
  vx = extract(data, "vx")
  vy = extract(data, "vy")
  trac = extract(data, "trac")
  x = extract(data, "x")*c_
  v2 = (vx**2+vy**2)
  lfac = 1./np.sqrt(1 - v2)
  ux = vx*lfac
  uy = vy*lfac
  dset = np.column_stack((rho, p, ux, uy, trac, x))  # cell values

  ntot = data.shape[0]
  Nr = pivot(data, "nact")[:,0] # number of cells in each track
  Nr = np.expand_dims(Nr,1)
  index = np.arange(ntot)
  data["index"] = index
  Index = pivot(data, "index")[:,0]
  Index = np.expand_dims(Index,1)
    # Index of first cell of track in memory
  th = pivot(data, "y")[:,0]
  dth = pivot(data, "dy")[:,0]
  t_jph = th - dth/2.
  t_jph = np.append(t_jph, t_jph[-1]+dth[-1]) # theta interface positions

  p_kph = np.array([0,2*np.pi]) # phi interface posisions

  return(dset, Index, Nr, T, p_kph, t_jph)



def toH5(key, it, index):

  if index > 9999:
    print("ERROR! index too high!")
    return -1

  pandata = readData(key, it, sequence=True)

  filename = '../../results/%s/%04d.h5'  %(key,index)
  if os.path.exists(filename):
    os.remove(filename)
  f = h5py.File(filename, 'a')

  values, Index, Nr, T, p_kph, t_jph = pandas2double(pandata)

  # fluid values
  data = f.create_group("Data")
  dset = f.create_dataset("Data/Cells", data=values)

  # Geometry
  grid = f.create_group("Grid")
  dset = f.create_dataset("Grid/Index", data=Index)
  dset = f.create_dataset("Grid/Nr", data=Nr)
  dset = f.create_dataset("Grid/T", data=T)
  dset = f.create_dataset("Grid/p_kph", data=p_kph)
  dset = f.create_dataset("Grid/t_jph", data=t_jph)

  f.close()

  return T



def prepForBlast(dir, geometry=None):
  if geometry is None or (geometry != "cylindrical" and geometry != "spherical"):
    print("ERROR! Please specify geometry (cylindrical/spherical).")
    print("Exiting and returning -1.")
    return(-1)
  header_filename = '../../results/%s/header.h5'  %(dir)
  if os.path.exists(header_filename):
    os.remove(header_filename)
  f = h5py.File(header_filename, 'a')

  utf8_type = h5py.string_dtype('utf-8', 12)
  dset_geom = f.create_dataset("geometry", 
                               data=np.array(geometry.encode("utf-8"), 
                                             dtype=utf8_type))
    # fixed-length string

  filelist, times = applyToAll(dir, toH5)
  no_snapshots = len(filelist)
  dset_no_snapshots = f.create_dataset("no_snapshots", data=no_snapshots)
  dset_t = f.create_dataset("t", data=times)
  f.close()














