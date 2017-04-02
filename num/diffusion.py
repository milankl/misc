## HEAT DIFFUSION

import numpy as np
import matplotlib.pyplot as plt
import pylan as pn

## grid

dx = 0.1
Nx = int(1/dx+1)
x = np.arange(0,1+dx,dx)

dt = 0.01
T = 5
Nt = int(T/dt)
t = np.arange(0,T,dt)

xx,tt = np.meshgrid(x,t)

## Laplace operator

a = (1/np.pi)**2 #alpha**2
L = a*(-2*np.eye(Nx) + np.diag(np.ones(Nx-1),1) + np.diag(np.ones(Nx-1),-1))/dx**2
L = L[1:-1,1:-1] #apply only in the interior

## Euler forward

def forward(T):
    Nt = int(T/dt)
    phi = np.zeros((Nt,Nx))
    phi[0,:] = np.cos(np.pi*(x-.5))

    for i in range(Nt-1):
        phi[i+1,1:-1] = phi[i,1:-1] + dt*L.dot(phi[i,1:-1])
    
    return phi

def backward(T):
    Nt = int(T/dt)
    phi = np.zeros((Nt,Nx))
    phi[0,:] = np.cos(np.pi*(x-.5))

    for i in range(Nt-1):
        phi[i+1,1:-1] = np.linalg.solve(np.eye(Nx-2) - dt*L,phi[i,1:-1])
    
    return phi

def analytic(T):
    return np.exp(-tt)*np.cos(np.pi*(xx-.5))

## plotting

ferror = forward(T)-analytic(T)
berror = backward(T)-analytic(T)

levmax = pn.round1(abs(np.hstack((ferror,berror))).max())*1.1

lev1 = np.linspace(0,1,11)
lev2 = np.linspace(-levmax,levmax,13)

fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)

im1 = ax1.contourf(xx,tt,forward(T),lev1)
ax2.contourf(xx,tt,backward(T),lev1)

im2 = ax3.contourf(xx,tt,ferror,lev2,cmap='RdBu_r')
ax4.contourf(xx,tt,berror,lev2,cmap='RdBu_r')


plt.colorbar(im1,ax=(ax1,ax2))
plt.colorbar(im2,ax=(ax3,ax4))

ax1.set_title('forward')
ax2.set_title('backward')
ax3.set_title('Error: forward - analytic')
ax4.set_title('Error: backward - analytic')

ax4.set_xlabel('x')
ax3.set_xlabel('x')
ax1.set_ylabel('t')
ax3.set_ylabel('t')

plt.show()