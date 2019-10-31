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


nodes = 13
T = 5           #total time
num_steps = 500 #number of steps
dt = T/num_steps #time step size


k = Constant(dt)
#set parameter values


# fluid viscosity
mu = 0.01

# stiffness parameter
K = 100
#reference length
l0 = 0.97 #micrometer 


    
def meantensionext(l):
    kt =0.012 #dyn/cm elastic modulus
    muext = 2*10**(-4) #dyn*s/cm viscosity
    return kt*(l/l0-1)+muext*1/l * ((l-ln)/k)

def meantensionint(L):
    muint = 10**(-4) #dyn*s/cm viscosity
    return muint*1/L * ((L-Ln)/k)

def meanshear(a,a1, l):
    #bending modulus on ith node
    kbending = 9*10**(-12) 
    return kbending*(a-an)/(l0*l)

# need to reevaluate alpha
    
def pint(A):
    #reference Area
    Aref = 22.2 #micrometer**2
    kp = 50 #dyn/cm pressure modulus
    return kp*(1-A/Aref)


X=np.array([])
Y=np.array([])

for i in range(nodes):
    alpha = 2/(nodes -1) *i *np.pi #alpha in rad
    x,y = 5, 5
    
    x_i = x +4 -4*np.cos(alpha)
    y_i = y +4* np.sin(alpha)
    
    X = np.append(X, x_i)
    Y = np.append(Y, y_i)
    
