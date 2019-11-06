# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 17:58:42 2019

@author: hanno
"""

import h5py
f = h5py.File('velocity2.h5', 'r')

print(f)
list(f.keys())