## NUMERICS 1

import numpy as np
import matplotlib.pyplot as plt

## constants
global H,f,taux,rho,nu
H = 100.
f = 1e-4
taux = 1e-1
rho = 1000.
nu = 1e-1

## analytic solution

beta = (1 + 1j)*np.sqrt(f/nu/2)
A = taux/nu/rho/beta*(1-1/(1-np.exp(2*beta*H)))
B = -taux/nu/rho/beta/(1-np.exp(2*beta*H))
za = np.linspace(-H,0,1000)
qa = (A*np.exp(beta*za) + B*np.exp(-beta*za))
ua = np.real(qa)
va = np.imag(qa)

##

def solve_ekman(N):
    
    #grid
    zb = np.linspace(0,H,N+1)
    dz = zb[1] - zb[0]
    z = np.arange(dz/2,H,dz)
    
    #laplace operator
    A = (np.diag(-np.ones(N)*2) +\
    np.diag(np.ones(N-1),-1) +\
    np.diag(np.ones(N-1),1))/dz**2
    
    #von Neumann boundary conditions
    A[0,0] = -1/dz**2
    A[-1,-1] = -1/dz**2
    
    #relaxation term
    B = -np.eye(N)*1j*f
        
    #source term
    s = np.zeros(N)
    s[0] = taux/rho/dz
    
    # solving
    q = np.linalg.solve(nu*A+B,-s)
    u = np.real(q)
    v = np.imag(q)
    return u,v,z

u10,v10,z10 = solve_ekman(10)
u5,v5,z5 = solve_ekman(5)
u20,v20,z20 = solve_ekman(20) 

## plotting

N = 10
zb = np.linspace(0,H,N+1)

fig,(ax1,ax2) = plt.subplots(1,2,sharey=True)
ax1.plot(u10,-z10,'ks-',label=r'$u$, dz=10m')

ax1.plot(ua,za,'r',label=r'$u_a$')
ax1.set_title(r'(a) zonal velocity $u$')
ax1.set_ylabel(r'depth $z$')
ax1.set_xlabel(r'[ms$^{-1}$]')

ax2.plot(v10,-z10,'ks-',label=r'$v$, dz=10m')
ax2.plot(va,za,'r',label=r'$v_a$')
ax2.set_title(r'(b) meridional velocity $v$')
ax2.set_xlabel(r'[ms$^{-1}$]')

ax1.plot(np.zeros(N),-z10,'+',color='grey',label='midpoints')
ax1.scatter(np.zeros(N+1),-zb,s=20,c='grey',marker=r'$-$',label='grid')
ax1.legend(loc=4)
ax1.set_ylim(-H,0)

ax2.plot(np.zeros(N),-z10,'+',color='grey',label='midpoints')
ax2.scatter(np.zeros(N+1),-zb,s=20,c='grey',marker=r'$-$',label='grid')
ax2.legend(loc=3)

fig2,(ax21,ax22) = plt.subplots(1,2,sharey=True)
ax21.plot(u5,-z5,'gs-',label=r'$u$, dz=20m')
ax21.plot(u20,-z20,'bs-',label=r'$u$, dz=5m')
ax21.plot(ua,za,'r',label=r'$u_a$')
ax21.set_title(r'(a) zonal velocity $u$')
ax21.set_ylabel(r'depth $z$')
ax21.set_xlabel(r'[ms$^{-1}$]')
ax21.legend(loc=4)

ax22.plot(v5,-z5,'gs-',label=r'$v$, dz=20m')
ax22.plot(v20,-z20,'bs-',label=r'$v$, dz=5m')
ax22.plot(va,za,'r',label=r'$v_a$')
ax22.set_title(r'(b) meridional velocity $v$')
ax22.set_xlabel(r'[ms$^{-1}$]')
ax22.legend(loc=3)

plt.show()