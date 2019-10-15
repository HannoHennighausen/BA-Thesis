#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 09:45:27 2019

@author: bq_hhennighausen
"""

from __future__ import print_function
from mshr import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


# Create mesh and define function spaces
mesh = Mesh()

# Create list of polygonal domain vertices
base = Rectangle (Point(0.0,0.0), Point(40.0, 10.0))
top = Rectangle (Point(10.0,6.5), Point(25.0, 10.0))
bottom = Rectangle (Point(10.0,0.0), Point(25.0, 3.5))

mesh = generate_mesh(base - top - bottom, 40)

plt.figure()
plot(mesh)
plt.show()