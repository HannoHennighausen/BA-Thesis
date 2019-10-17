#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:39:33 2019

@author: bq_hhennighausen
"""

from __future__ import print_function
from mshr import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

mesh = Mesh()


base = Rectangle (Point(0.0,0.0), Point(40.0, 10.0))
top = Rectangle (Point(10.0,6.5), Point(25.0, 10.0))
bottom = Rectangle (Point(10.0,0.0), Point(25.0, 3.5))
cylinder = Circle(Point(15.0, 5.0), 0.5)

mesh = generate_mesh(base - top - bottom - cylinder, 400)



# Define function spaces
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

# Define boundaries
inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 40)'
def wall_1(x, on_boundary):
    return on_boundary and (between(x[0], (0.0,10.0)) and (near(x[1], 0) or near(x[1], 10)))
def wall_2(x, on_boundary):
    return on_boundary and (between(x[0], (10.0,25.0)) and (near(x[1], 6.5) or near(x[1], 3.5)))
def wall_3(x, on_boundary):
    return on_boundary and (between(x[0], (25.0,40.0)) and (near(x[1], 0) or near(x[1], 10)))
def wall_4(x, on_boundary):
    return on_boundary and (between(x[1], (0.0,3.5)) and (near(x[0], 10) or near(x[0], 25)) )
def wall_5(x, on_boundary):
    return on_boundary and (between(x[1], (6.5,10.0)) and (near(x[0], 10) or near(x[0], 25)))
cylinder = 'on_boundary && x[0]>14.4 && x[0]<15.6 && x[1]>4.4 && x[1]<5.6'

# Define inflow profile
#inflow_profile = ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0')

# Define boundary conditions
bcu_inflow = DirichletBC(Q, Constant(8), inflow)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu_noslip1  = DirichletBC(V, Constant((0, 0)), wall_1)
bcu_noslip2  = DirichletBC(V, Constant((0, 0)), wall_2)
bcu_noslip3  = DirichletBC(V, Constant((0, 0)), wall_3)
bcu_noslip4  = DirichletBC(V, Constant((0, 0)), wall_4)
bcu_noslip5  = DirichletBC(V, Constant((0, 0)), wall_5)
bcu = [ bcu_noslip1, bcu_noslip2, bcu_noslip3, bcu_noslip4, bcu_noslip5, bcu_cylinder]
bcp = [bcu_inflow, bcp_outflow]

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)


f = Constant(( 0.0, 0.0))
a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx
L = inner(f, v)*dx

# Form for use in constructing preconditioner matrix
b = inner(grad(u), grad(v))*dx + p*q*dx

# Assemble system
A, bb = assemble_system(a, L, bcs)

# Assemble preconditioner system
P, btmp = assemble_system(b, L, bcs)

# Create Krylov solver and AMG preconditioner
solver = KrylovSolver(krylov_method, "amg")

# Associate operator (A) and preconditioner matrix (P)
solver.set_operators(A, P)

# Solve
U = Function(W)
solver.solve(U.vector(), bb)

# Get sub-functions
u, p = U.split()


# Plot solution
plt.figure()
plot(u)
plot(p)
interactive()
plt.show()