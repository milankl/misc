## NUMERICS ASSIGNMENT 4
# ADVECTION
""" Solving du/dt + du/dx = 0 """

import numpy as np
import matplotlib.pyplot as plt
import pylan as pn

## set up the grid
dx = .5
xend = 30
x = np.arange(0,xend+dx,dx)
nx = len(x)

tend = 10
## analytical solution
def ua(t):
    return (((x-t) >= 1.5) * ((x-t) <= 6.5))*1

## set up LEFT-DERIVATIVE, CENTRED DERIVATIVE
Gl = (np.eye(nx) + np.diag(-np.ones(nx-1),-1))/dx
Gc = (np.diag(np.ones(nx-1),1) + np.diag(-np.ones(nx-1),-1))/2/dx

#Gl[0,-1] = -1
#Gc[0,-1] = -1
#Gc[-1,0] = 1

## forward in time, centred in space
def FTUP(nt,dt):
    u = np.zeros((nt,nx))
    u[0,:] = ((x >= 1.5) * (x <= 6.5))*1 #initial conditions
    
    for i in range(nt-1):
        u[i+1,:] = u[i,:] - dt*Gl.dot(u[i,:])
    
    return u

## leap frog in time, centered in space    
def LFCS(nt,dt):
    u = np.zeros((nt,nx))
    u[:2,:] = FTUP(2,dt) #use forward in time for n = 1
    
    for i in range(1,nt-1):
        u[i+1,:] = u[i-1,:] - 2*dt*Gc.dot(u[i,:])
    
    return u

## plotting
dtlist = np.array([.1,.5,1])

fig,axs = plt.subplots(2,3,sharex=True,sharey=True,figsize=(11,6))
for i in range(3):
    
    n = int(tend/dtlist[i])+1
    u = FTUP(n,dtlist[i])
    u2 = LFCS(n,dtlist[i])
    
    axs[0,i].plot(x,u[int(n/2),:],label='FT+UP')
    axs[0,i].plot(x,u2[int(n/2),:],label='LF+CS')
    axs[0,i].plot(x,ua(5.),'--',label='analytic',lw=2)
    
    axs[1,i].plot(x,u[-1,:])
    axs[1,i].plot(x,u2[-1,:])
    axs[1,i].plot(x,ua(10.),'--',lw=2)

    if i == 0:
        axs[0,i].legend()
        axs[0,i].set_ylabel(r'$u$')
        axs[1,i].set_ylabel(r'$u$')
    
    axs[1,i].set_xlabel(r'$x$')
    axs[0,i].set_title(r'$t = %i, \Delta t = %.1f$' % (int(tend/2),dtlist[i]))
    axs[1,i].set_title(r'$t = %i, \Delta t = %.1f$' % (int(tend),dtlist[i]))
    

axs[0,0].set_ylim(-.5,1.5)

plt.tight_layout()
plt.show()