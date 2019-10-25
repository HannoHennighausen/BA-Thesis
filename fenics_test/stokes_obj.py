#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:50:38 2019

@author: bq_hhennighausen

FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for flow around a cylinder using the Incremental Pressure Correction
Scheme (IPCS).
  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

from __future__ import print_function
from mshr import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

T = 3.0         # final time
num_steps = 100   # number of time steps
dt = T / num_steps # time step size
mu = 0.001       # dynamic viscosity
rho = 1            # density
e = 0.01        #epsilon
# Create mesh
# Create mesh and define function spaces
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

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx \
   + e * rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Create XDMF files for visualization output
xdmffile_u = XDMFFile('cylinder/velocity2.xdmf')
xdmffile_p = XDMFFile('cylinder/pressure2.xdmf')

# Create time series (for use in reaction_system.py)
timeseries_u = TimeSeries('cylinder/velocity_series2')
timeseries_p = TimeSeries('cylinder/pressure_series2')

# Save mesh to file (for use in reaction_system.py)
File('cylinder/cylinder2.xml.gz') << mesh

# Create progress bar
#progress = Progress('Time-stepping')
#set_log_level(PROGRESS)

# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg', 'sor')

    # Plot solution
    plot(u_, title='Velocity')
    plot(p_, title='Pressure')

    

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

    # Update progress bar
    #progress.update(t / T)
    print('u max:', u_.vector().get_local().max())


# Save solution to file (XDMF/HDF5)
xdmffile_u.write(u_, t)
xdmffile_p.write(p_, t)

# Save nodal values to file
timeseries_u.store(u_.vector(), t)
timeseries_p.store(p_.vector(), t)


# Hold plot
#interactive()

