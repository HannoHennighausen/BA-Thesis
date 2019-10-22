#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:30:17 2019

@author: bq_hhennighausen
"""

from __future__ import print_function
from mshr import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


from dolfin import *

# Load mesh and subdomains

# Create mesh and define function spaces
mesh = UnitSquareMesh(8 , 8)
'''
# Create list of polygonal domain vertices
base = Rectangle (Point(0.0,0.0), Point(40.0, 10.0))
top = Rectangle (Point(10.0,6.5), Point(25.0, 10.0))
bottom = Rectangle (Point(10.0,0.0), Point(25.0, 3.5))
cylinder  = Circle(Point(5.0, 5.0), 4.0)
#cylinder  = Ellipse(Point(18.0,5.0), 6.0,1.5)
#cylinder  = Ellipse(Point(32.0,5.0), 5.0,3.0)
domain = base - top - bottom - cylinder
#domain.set_subdomain(1, cylinder)
mesh = generate_mesh(domain, 50);
'''
# Define function spaces
V = VectorElement('Lagrange',mesh.ufl_cell(), 2)
Q = FiniteElement('Lagrange',mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, V*Q)

# Define boundaries
inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 8)'
def wall_1(x, on_boundary):
    return on_boundary and (near(x[1], 0) or near(x[1], 8))



# Define boundary conditions
bcp_inflow = DirichletBC(W.sub(1), Constant(1), inflow)
bcp_outflow = DirichletBC(W.sub(1), Constant(0), outflow)
bcu_noslip1  = DirichletBC(W.sub(0), Constant((0, 0)), wall_1)


#bcu = [ bcu_noslip1, bcu_noslip2, bcu_noslip3, bcu_noslip4, bcu_noslip5, bcu_cylinder]
#bcp = [bcp_inflow, bcp_outflow]
bcs = [bcu_noslip1, bcp_inflow, bcp_outflow]

# Define trial and test functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

f = Constant((0, 0))
a = inner(grad(u), grad(v))*dx + p*div(v)*dx + q*div(u)*dx
L = inner(f, v)*dx
# Compute solution
w = Function(W)
solve(a == L, w, bcs)
u, p = w.split()




# Save solution in VTK format
ufile_pvd = File("velocity9.pvd")
ufile_pvd << u
pfile_pvd = File("pressure9.pvd")
pfile_pvd << p

# Plot solution
plt.figure()
plot(u)
plt.show()

plt.figure()
plot(p)
plt.show()