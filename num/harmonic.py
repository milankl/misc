## NUMERICS HARMONIC DIFFUSION

import numpy as np
import matplotlib.pyplot as plt

## constants
A2,A4 = 1.,1.
dx = 1.
kmax = np.pi/dx
k = np.linspace(1e-2,kmax,1000)

tau2 = -dx**2/(2*A2*(np.cos(k*dx)-1))
tau4 = -dx**4/(2*A4*(4*np.cos(k*dx) - np.cos(2*k*dx)-3))

tau2a = 1./A2/k**2
tau4a = 1./A4/k**4

## plotting

fig,ax1 = plt.subplots()

ax1.plot(k,tau2,'g--',label=r'$\tau_{2,discr.}$',lw=2)
ax1.plot(k,tau2a,'g-',label=r'$\tau_{2,analytic}$',lw=2)

ax1.plot(k,tau4,'r--',label=r'$\tau_{4,discr.}$',lw=2)
ax1.plot(k,tau4a,'r-',label=r'$\tau_{4,analytic}$',lw=2)

ax1.set_yscale('log')
ax1.plot([kmax, kmax],[1e-2,1e3],'grey',ls='--',label=r'$k_{max}$')
ax1.set_xlabel(r'$|k|$')
ax1.set_ylabel(r'Damping time scale $\tau$')
ax1.set_ylim(1e-2,1e3)
ax1.set_title('Discrete (bi-)harmonic diffusion: damping time scales')
ax1.legend(loc=9)


plt.tight_layout()
plt.show()