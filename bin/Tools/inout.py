# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2021-01-30 17:59:32
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-01-31 16:12:43

import pandas as pd
import numpy as np

def readData(key, it=None, sequence=False):
  if sequence:
    filename = '../../results/%s/phys%010d.out'  %(key,it)
  elif it==None:
    filename = '../../results/%s'  %(key)
  else:
    filename = '../../results/%s%d.out'  %(key,it)
  data = pd.read_csv(filename, sep=" ")
  return(data)

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
