## WAVE EQUATION

import numpy as np
import matplotlib.pyplot as plt
import pylan as pn

## constants and grid
H = 10
L = 1e5
g = 9.8
F = 0.01/1e3/H #tau/rho0/H

dx = 5e3
dt = 300

cfl = np.sqrt(g*H)*dt/dx
print('cfl = %1.3f' % cfl)

T = 48*3600
N = int(T/dt)+1

## staggered grid
xu = np.arange(0,L+dx,dx)
xe = xu[:-1]+dx/2

nxu = xu.shape[0] - 2
nxe = xe.shape[0]

t = np.arange(0,T+dt,dt)
xxu,ttu = np.meshgrid(xu,t)
xxe,tte = np.meshgrid(xe,t)

## dx gradients
Gxu = (np.diag(np.ones(nxu+1),0) - np.diag(np.ones(nxe-1),-1))[:,:-1]/dx
Gxe = -Gxu.T

## preallocate

u = np.zeros((N,nxu))
eta = np.zeros((N,nxe))

for i in range(N-1):
    eta[i+1,:] = eta[i,:] - dt*H*Gxu.dot(u[i,:])
    u[i+1,:] = u[i,:] - g*dt*Gxe.dot(eta[i+1,:]) + dt*F

#pad u with zeros
u = np.concatenate((np.zeros((N,1)),u,np.zeros((N,1))),axis=1)

## plotting

fig,(ax1,ax2) = plt.subplots(1,2,sharex=True,sharey=True)

c1 = ax1.contourf(xxu/1e3,ttu/3600,u*1e2)
ax1.set_xlabel(r'$x$ [km]')
ax1.set_ylabel(r'$t$ [h]')
ax1.set_title(r'$u$ [cm/s]')
plt.colorbar(c1,ax=ax1)

c2 = ax2.contourf(xxe/1e3,tte/3600,eta*1e2)
ax2.set_xlabel(r'$x$ [km]')
ax2.set_title(r'$\eta$ [cm]')
plt.colorbar(c2,ax=ax2)

plt.show()




