# -*- coding: utf-8 -*-
# @Author: Eliot Ayache
# @Date:   2021-01-30 17:59:32
# @Last Modified by:   Eliot Ayache
# @Last Modified time: 2021-01-30 18:17:21

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
