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


from dolfin import *


# Test for PETSc or Tpetra
if not has_linear_algebra_backend("PETSc") and not has_linear_algebra_backend("Tpetra"):
    info("DOLFIN has not been configured with Trilinos or PETSc. Exiting.")
    exit()

if not has_krylov_solver_preconditioner("amg"):
    info("Sorry, this demo is only available when DOLFIN is compiled with AMG "
	 "preconditioner, Hypre or ML.")
    exit()

if has_krylov_solver_method("minres"):
    krylov_method = "minres"
elif has_krylov_solver_method("tfqmr"):
    krylov_method = "tfqmr"
else:
    info("Default linear algebra backend was not compiled with MINRES or TFQMR "
         "Krylov subspace method. Terminating.")
    exit()

# Load mesh and subdomains

# Create mesh and define function spaces
mesh = Mesh()

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

# Define function spaces

V = VectorElement('P',mesh.ufl_cell(), 1)
Q = FiniteElement('P',mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, V*Q)

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
cylinder = 'on_boundary && x[0]>0.9 && x[0]<9.1 && x[1]>0.9 && x[1]<9.1'

# Define inflow profile
#inflow_profile = ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0')

# Define boundary conditions
bcp_inflow = DirichletBC(W.sub(1), Constant(15), inflow)
bcp_outflow = DirichletBC(W.sub(1), Constant(0), outflow)
bcu_cylinder = DirichletBC(W.sub(0), Constant((0, 0)), cylinder)
bcu_noslip1  = DirichletBC(W.sub(0), Constant((0, 0)), wall_1)
bcu_noslip2  = DirichletBC(W.sub(0), Constant((0, 0)), wall_2)
bcu_noslip3  = DirichletBC(W.sub(0), Constant((0, 0)), wall_3)
bcu_noslip4  = DirichletBC(W.sub(0), Constant((0, 0)), wall_4)
bcu_noslip5  = DirichletBC(W.sub(0), Constant((0, 0)), wall_5)
#bcu = [ bcu_noslip1, bcu_noslip2, bcu_noslip3, bcu_noslip4, bcu_noslip5, bcu_cylinder]
#bcp = [bcp_inflow, bcp_outflow]
bcs = [bcu_noslip1, bcu_noslip2, bcu_noslip3, bcu_noslip4, bcu_noslip5, bcu_cylinder, bcp_inflow, bcp_outflow]

# Define trial and test functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

f = Constant((0, 0))
a = 1*(inner(grad(u), grad(v))*dx) + div(v)*p*dx + q*div(u)*dx
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
print(A)
# Solve
U = Function(W)
solver.solve(U.vector(), bb)

# Get sub-functions
(u, p) = U.split()



# Save solution in VTK format
ufile_pvd = File("velocity123.pvd")
ufile_pvd << u
pfile_pvd = File("pressure123.pvd")
pfile_pvd << p

# Plot solution
plt.figure()
plot(u)
plt.show()

plt.figure()
plot(p)
plt.show()