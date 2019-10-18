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
#cylinder  = Circle(Point(5.0, 5.0), 4.0)
#cylinder  = Ellipse(Point(18.0,5.0), 6.0,1.5)
cylinder  = Ellipse(Point(32.0,5.0), 5.0,3.0)
domain = base - top - bottom - cylinder
domain.set_subdomain(1, cylinder)
mesh = generate_mesh(domain, 50)

plt.figure()
plot(mesh)
plt.savefig('mesh3.pdf')
plt.show()
