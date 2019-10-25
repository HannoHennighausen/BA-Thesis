#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:24:46 2019

@author: bq_hhennighausen
"""

from __future__ import print_function
from mshr import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#set parameter values
# fluid viscosity
mu = 0.01

# stiffness parameter
K = 100
#reference length
l_0 = 0.97 #micrometer 

#external elements
k_t =0.012 #dyn/cm elastic modulus
mu_e = 2*10**(-4) #dyn*s/cm viscosity

#bending modulus on ith node
k_b = 9*10**(-12) 

#internal elements
mu_i = 10**(-4)

#reference Area
A_ref = 22.2 #micrometer**2
k_p = 50 #dyn/cm pressure modulus

