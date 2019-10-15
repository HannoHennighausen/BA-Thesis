#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 11:37:36 2019

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


#domain vertices
domain_vertices=[Point(0.0,0.0), 
                 Point(10.0,0.0),
                 Point(10.0,3.5),
                 Point(25.0,3.5),
                 Point(25.0,0.0), 
                 Point(40.0,0.0),
                 Point(40.0,10.0),
                 Point(25.0,10.0),
                 Point(25.0,6.5), 
                 Point(10.0,6.5),
                 Point(10.0,10.0), 
                 Point(0.0,10.0),
                 Point(0.0,0.0) ]
domain = Polygon(domain_vertices)
mesh = generate_mesh(domain,40)

# Save mesh to file (for use in reaction_system.py)
xdmffile_mesh = XDMFFile('Mesh/mesh.xdmf')
File('Mesh/mesh.xml.gz') << mesh

plt.figure()
plot(mesh)
plt.show()
'''
# Define boundaries
inflow  = 'near(x[0], 0)'
outflow = 'near(x[0], 1)'
walls   = 'near(x[1], 0) || near(x[1], 1)'

# Define boundary conditions
bcu_noslip  = DirichletBC(V, Constant((0, 0)), walls)
bcp_inflow  = DirichletBC(Q, Constant(8), inflow)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_noslip]
bcp = [bcp_inflow, bcp_outflow]
'''
