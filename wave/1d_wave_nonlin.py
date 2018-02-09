## WAVE EQ
""" this script solves the 1D-wave equations, i.e.
    
        du/dt = -d/dx(1/2*u^2 + g*h - nu*du/dx) + F
        dh/dt = -d/dx(u*h)."""
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time as tictoc

## physical parameters
Lx = 1.     # D = (0,Lx) is the domain
H = 10.      # depth at rest
g = 10.     # gravitational acceleration

## numerical parameters
nx = 500                    # number of grid points for the variable h
nu = nx                   # number of grid points for u
nt = 3000                     # number of time steps to integrate

cfl = 0.9                   # cfl number, 1 is the limit

dx = Lx/nx                  # grid spacing

# staggered grid
x_h = np.arange(dx/2.,Lx,dx)    # grid x on h-points
x_u = np.arange(dx,Lx-dx+dx/2.,dx)    # on u-points

cph = np.sqrt(g*H)          # long wave phase speed
dt = cfl * dx / cph           # dt based on desired cfl-number

# parameters
rho = 1
v = 0.006       # dependent on the grid spacing! choose carefully
F0 = 500.
F = -F0*np.sin(x_h/Lx*8*np.pi)*np.sin(x_h/Lx*np.pi)**2

## OPERATORS
# x gradient on the T-grid
GTx = (sparse.diags(np.ones(nx),0,shape=(nx,nx)) +\
    sparse.diags(-np.ones(nx),-1,shape=(nx,nx))) / dx

GTx[0,-1] = -1/dx   # for periodic boundary conditions

# x gradient on the u-grid
Gux = -GTx.T

# interpolation from u to T and vice versa
ITu = abs(GTx*dx/2.)
IuT = abs(Gux*dx/2.)

## initial conditions
h0 = np.zeros(nx)
u0 = np.zeros(nu)

##  time integration - TRIVIAL WAY
h1 = h0.copy()        # store initial conditions
u1 = u0.copy() 

def rhs(u,h):
    # h is actually eta, h_u is eta+H on the u-grid
    h_u = ITu.dot(h+H)
    
    # non-linear
    rhs_u = -GTx.dot(.5*IuT.dot(u**2) + g*h - v*Gux.dot(u)) + 1/rho/h_u*F0*np.sin(x_h/Lx*8*np.pi - 0.05*it)*np.sin(x_h/Lx*np.pi)**2
    rhs_h = -Gux.dot(u*h_u)
    
    return rhs_u,rhs_h

def RK4(u,h):
    k1 = rhs(u,h)
    k2 = rhs(u + dt/2.*k1[0],h + dt/2.*k1[1])
    k3 = rhs(u + dt/2.*k2[0],h + dt/2.*k2[1])
    k4 = rhs(u + dt*k3[0],h + dt*k3[1])
    
    du = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.
    dh = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.
    return du,dh

fig,(ax1,ax2,ax3) = plt.subplots(3,1,figsize=(9,8))

mass1 = np.empty(nt)
energy1 = np.empty(nt)

time = np.arange(nt)*dt

for it in range(nt):
    
    mass1[it] = (h1+H).sum()/(H*nx)
    energy1[it] = (.5*(h1+H)*IuT.dot(u1**2) + .5*h1**2).mean()
    
    if it == 0:
        l11, = ax1.plot(x_h,h1)
        ax1.set_ylim(H-2,H+2)
        ax1.set_xlabel('x')
        ax1.set_ylabel(r'$h$')
        ax1.set_title(r'surface displacement $h$')
        
        l21, = ax2.plot(time[:it],mass1[:it])
        ax2.set_ylim(0.95,1.05)
        ax2.plot([0,nt*dt],[mass1[0],mass1[0]],'grey',lw=0.5)
        ax2.set_xlabel('time')
        ax2.set_ylabel('mass')
        
        l31, = ax3.plot(time[:it],energy1[:it])
        ax3.set_ylim(0,4)
        ax3.set_xlim(ax2.get_xlim())
        ax3.set_xlabel('time')
        ax3.set_ylabel('energy')
        
        plt.tight_layout()
        
    else:
        du1,dh1 = RK4(u1,h1)
        u1 += dt*du1
        h1 += dt*dh1
        
        if it % 5 == 1:
            plt.pause(0.00001)
            l11.set_data(x_h,h1+H)
            l21.set_data(time[:it],mass1[:it])
            l31.set_data(time[:it],energy1[:it])
            
plt.pause(2)
plt.close(fig)
        
        


