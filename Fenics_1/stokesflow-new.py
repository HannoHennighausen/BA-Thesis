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

# Create mesh and define function spaces
mesh = Mesh()

# Create microchannel as mesh
base = Rectangle (Point(0.0,0.0), Point(40.0, 10.0))
top = Rectangle (Point(10.0,6.5), Point(25.0, 10.0))
bottom = Rectangle (Point(10.0,0.0), Point(25.0, 3.5))
cylinder  = Circle(Point(5.0, 5.0), 4.0)
#cylinder  = Ellipse(Point(18.0,5.0), 6.0,1.5)
cylinder1  = Ellipse(Point(32.0,5.0), 5.0,3.0)
domain = base - top - bottom - cylinder - cylinder1

mesh = generate_mesh(domain, 100);

# Define function spaces

#2dim for velocity
V = VectorElement('P',mesh.ufl_cell(), 2)
VV= FunctionSpace(mesh, V)

#1dim for pressure
Q = FiniteElement('P',mesh.ufl_cell(), 1)
QQ= FunctionSpace(mesh, Q)

#create vectorspace by multipling 2x1 dimensions
W = FunctionSpace(mesh, V*Q)

# Define boundaries


def left(x, on_boundary): return x[0] < DOLFIN_EPS
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
#cylinder = 'on_boundary && x[0]>27 && x[0]<37 && x[1]>2 && x[1]<8'
cylinder1 = 'on_boundary && x[0]>26.9 && x[0]<37.1 && x[1]>2 && x[1]<8'



# Define boundary conditions

# Define inflow profile
inflow = Expression(('(1-(x[1]-5)*(x[1]-5)/25)', '0.0'), degree=2)
outflow  = 'near(x[0], 40)'
bcp_inflow = DirichletBC(W.sub(0),  inflow, left)
bcp_outflow = DirichletBC(W.sub(1), Constant(0), outflow)
bcu_cylinder = DirichletBC(W.sub(0), Constant((0, 0)), cylinder)
bcu_cylinder1 = DirichletBC(W.sub(0), Constant((0, 0)), cylinder1)
bcu_noslip1  = DirichletBC(W.sub(0), Constant((0, 0)), wall_1)
bcu_noslip2  = DirichletBC(W.sub(0), Constant((0, 0)), wall_2)
bcu_noslip3  = DirichletBC(W.sub(0), Constant((0, 0)), wall_3)
bcu_noslip4  = DirichletBC(W.sub(0), Constant((0, 0)), wall_4)
bcu_noslip5  = DirichletBC(W.sub(0), Constant((0, 0)), wall_5)

# collect boundary conditions

bcs = [bcu_noslip1, bcu_noslip2, bcu_noslip3, bcu_noslip4, bcu_noslip5, bcu_cylinder, bcp_inflow,\
       bcp_outflow, bcu_cylinder1]

# Define trial and test functions
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)


#define variational problem for Stokes flow
f = Constant((0, 0))
a = 0.1*(inner(grad(u), grad(v))*dx) + div(v)*p*dx + q*div(u)*dx
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
(u, p) = U.split()


#extract velocity values
u_P1 = Function(VV)
u_P1 = project(u, VV)
u_nodal_values = u_P1.vector()
u_array = u_nodal_values.vec()

#print(u_array.array)
#print(u_array[0])
#
i=2
coor = mesh.coordinates()
#print(coor.shape[0])
#print(coor[i][0], coor[i][1], u_array[i])

#extract pressure values
p_P1 = Function(QQ)
p_P1 = project(p, QQ)
p_nodal_values = p_P1.vector()
p_array = p_nodal_values.vec()

# print(p_array.array)
# print(p_array.size)
vertex_values = u_P1[0].compute_vertex_values()
for i, x in enumerate(coor):
    print('vertex %d: vertex_values[%d] = %g\tu(%s) = %g' %(i, i, vertex_values[i], x, u(x)))

arr = u.vector().get_local()
coor = mesh.coordinates()
vtd = W.dofmap()

values = list()
#print(coor)
np.savetxt('coor.txt', coor, fmt='%s')
np.savetxt('velocities.txt',u_array)
#print(str(vtd.entity_dofs(mesh, 0)))
# for i, dum in enumerate(coor):
#     values.append([arr[vtd[2*i]],arr[vtd[2*i+1]]])
# values = np.array(values)
#
# x = list()
# y = list()
# vx = list()
# vy = list()
# for i, dum in enumerate(coor):
#     print ('(%f,%f) -> (%f,%f)' %(coor[i][0], coor[i][1], values[i][0], values[i][1]))
#     x.append(coor[i][0])
#     y.append(coor[i][1])
#     vx.append(values[i][0])
#     vy.append(values[i][1])
#
# print(x,y, vx,vy)
# Save solution in VTK format
ufile_pvd = File("velocity_front and back.pvd")
ufile_pvd << u
pfile_pvd = File("pressure_front and back.pvd")
pfile_pvd << p

# Plot solution
#plt.figure()
#plot(u)
#plt.show()

#plt.figure()
#plot(p)
#plt.show()

Wortis