## WAVE EQ
""" this script solves the 1D-wave equations, i.e.
    
        du/dt = -g*dh/dx
        dh/dt = -H*du/xu."""
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time as tictoc

## physical parameters
Lx = 1.     # D = (0,Lx) is the domain
H = 1.      # depth at rest
g = 10.     # gravitational acceleration

## numerical parameters
nx = 500000                    # number of grid points for the variable h
nu = nx-1                   # number of grid points for u
nt = 10                     # number of time steps to integrate

cfl = 0.9                   # cfl number, 1 is the limit

dx = Lx/nx                  # grid spacing

x_h = np.arange(dx/2.,Lx,dx)    # grid x on h-points
x_u = np.arange(dx,Lx-dx,dx)    # on u-points

cph = np.sqrt(g*H)          # long wave phase speed
dt = cfl * dx / cph           # dt based on desired cfl-number

## OPERATORS
GTx = (sparse.diags(np.ones(nx-1),1,shape=(nu,nx)) +\
    sparse.diags(-np.ones(nx-1),0,shape=(nu,nx))) / dx
    
Gux = -GTx.T

## initial conditions
# gaussian initial state
sigma = 0.08
h0 = np.exp(-(x_h-.4)**2/(2*sigma**2))
u0 = np.zeros(nu)

##  time integration - TRIVIAL WAY
h = h0.copy()        # store initial conditions
u = u0.copy() 

def rhs(u,h):
    dhdx = GTx.dot(h)
    dudx = Gux.dot(u)
    
    rhs_u = -g*dhdx
    rhs_h = -H*dudx
    
    return rhs_u,rhs_h

def RK4(u,h):
    k1 = rhs(u,h)
    k2 = rhs(u + dt/2.*k1[0],h + dt/2.*k1[1])
    k3 = rhs(u + dt/2.*k2[0],h + dt/2.*k2[1])
    k4 = rhs(u + dt*k3[0],h + dt*k3[1])
    
    du = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.
    dh = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.
    return du,dh

tic = tictoc.time()
for i in range(nt):
        du,dh = RK4(u,h)
        u += dt*du
        h += dt*dh

t1 = tictoc.time() - tic

##  time integration - OPTIMIZED WAY
h = h0.copy()              # store initial conditions
u = u0.copy() 

# preallocate rhs
du = np.empty_like(u)
dh = np.empty_like(h)

# preallocate RK4 memories
u1 = u.copy()
u2 = u.copy()

h1 = h.copy()
h2 = h.copy()

a = np.array([1./6.,1./3.,1./3.,1./6.])
b = np.array([0.5,0.5,1.])

def rhs2(du,dh,dhdx,dudx):
    #dhdx[:] = 
    #dudx[:] = 
    
    #du[:] = 
    #dh[:] = 
    return du,dh

tic = tictoc.time()
for i in range(nt):
        
        u1[:] = u
        h1[:] = h
        
        for rk in range(4):  
            du = -g*GTx.dot(h)
            dh = -H*Gux.dot(u)
        
            if rk < 3:
                u = u1 + b[rk]*dt*du
                h = h1 + b[rk]*dt*dh
        
            u2 += a[rk]*dt*du
            h2 += a[rk]*dt*dh
        
        u[:] = u2
        h[:] = h2

t2 = tictoc.time() - tic
print(t1/t2)
"""
# for comparison
u_opt = u.copy()
h_opt = h.copy()

plt.plot(h0)
plt.plot(h)

plt.show()"""