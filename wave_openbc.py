## WAVE EQ
""" this script solves the 1D-wave equations, i.e.
    
        du/dt = -g*dh/dx
        dh/dt = -H*du/xu
    
    in the following form with h as prognostic variable
    
        d/dt^2 h = g*H*(d/dx)^2 h
        
    with a Laplacian in space and also a Laplacian in time.
    The boundary conditions are Dirichlet forced in the left and open in the right."""
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

## physical parameters
Lx = 1.     # D = (0,Lx) is the domain
H = 1.      # depth at rest
g = 10.     # gravitational acceleration

## numerical parameters
nx = 400                    # number of grid points for the variable h
nt = 800                    # number of time steps to integrate

cfl = 1.0                   # cfl number, 1 is the limit

dx = Lx/nx                  # grid spacing
x = np.arange(dx/2,Lx,dx)   # grid x
cph = np.sqrt(g*H)          # long wave phase speed
dt = cfl * dx / cph           # dt based on desired cfl-number
q = (1-cfl)/(1+cfl)         # for open boundary conditions

## Laplace operator in space - 2nd order in space the (1,-2,1)-stencil
# with dirichlet boundary conditions, without the 1/dx**2 factor
L = (sparse.diags(-2*np.ones(nx),0,shape=(nx,nx)) +\
    sparse.diags(np.ones(nx),1,shape=(nx,nx)) +\
    sparse.diags(np.ones(nx),-1,shape=(nx,nx))).tocsr()

# von neumann conditions at the right boundary for open boundary conditions
L[-1,-1] = -1

## boundary conditions at the right boundary - as source term s
omega = 40
def s(t):
    z = np.zeros(nx)
    z[0] = np.sin(omega*t)
    return z

## initial conditions
# gaussian initial state
#sigma = 0.08
#h0 = np.exp(-(x-.4)**2/(2*sigma**2))

# rest
h0 = np.zeros_like(x)

##  time integration - Laplacian in time
h = np.zeros((nt,nx))    # preallocate
h[0,:] = h0              # store initial conditions
t = 0

# the first time step, use second intial condition: d/dt^2 h = 0 for t = 0
h[1,:] = h[0,:] + cfl**2*(L.dot(h[0,:]) + s(t))
t += dt

for i in range(1,nt-1):
    # Laplacian in time
    h[i+1,:] = 2*h[i,:] - h[i-1,:] + cfl**2*(L.dot(h[i,:]) + s(t))
    
    # predict right-most point based on advection equation to yield open boundary conditions
    h[i+1,-1] = h[i,-2] - q*h[i+1,-2] + q*h[i,-1]
    
    # update time
    t += dt

## PLOTTING
ax = plt.subplots(1,1,figsize=(9,4))

for frame in range(nt):
    if frame == 0:
        l1, = plt.plot(x,h[frame,:])
        plt.ylim(-1.5,1.5)
    else:
        plt.pause(0.00001)
        l1.set_data(x,h[frame,:])

plt.pause(1)
plt.close('all')