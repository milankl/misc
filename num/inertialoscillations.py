## INERTIAL OSCILLATIONS

import numpy as np
import matplotlib.pyplot as plt

## constants

f = 1e-4
q0 = (1+1j)/np.sqrt(2)
T = 10*2*np.pi/f

## functions
def forward(dt): #euler forward
    return (1-1j*f*dt)**np.arange(int(T/dt))*q0

def backward(dt):
    return (1/(1+1j*f*dt))**np.arange(int(T/dt))*q0

def leapfrog(dt):
    N = int(T/dt)
    q = np.zeros(N,dtype=complex)
    q[0] = q0
    q[1] = backward(dt)[1] #use backward for first step
    
    for i in range(1,N-1):
        q[i+1] = q[i-1] - 2*1j*f*q[i]*dt
    
    return q
    
def leapfrog2(dt):
    n = np.arange(int(T/dt))
    lam1 = np.sqrt(1-(f*dt)**2) - 1j*f*dt
    return lam1**n*q0
    
def ana(dt): #analytic
    t = np.arange(0,T/10,dt/10)
    return q0*np.exp(-1j*f*t)

## evaluate

dtlist = [100,200,500,1000]

fig,ax = plt.subplots(2,2,sharex=True,sharey=True)

for i in range(len(dtlist)):
    
    qf = forward(dtlist[i])
    qb = backward(dtlist[i])
    ql = leapfrog(dtlist[i])
    qa = ana(dtlist[i])
    
    j = int(i/2) % 2
    k = i % 2
    
    ax[j,k].plot(np.real(qa),np.imag(qa),'r',lw=5,label='analytic',alpha=.5)
    ax[j,k].plot(np.real(qb),np.imag(qb),'b',label='backward')
    ax[j,k].plot(np.real(qf),np.imag(qf),'g',label='forward')
    ax[j,k].plot(np.real(ql),np.imag(ql),'k',label='leapfrog')
    ax[j,k].set_title(r'$\Delta t$ = %is' % dtlist[i])
    ax[j,k].scatter(np.real(q0),np.imag(q0),50,'#FF0000',marker='v',zorder=5)
    ax[j,k].scatter(np.real(qb[-1]),np.imag(qb[-1]),50,'b',marker='v',zorder=5)
    ax[j,k].scatter(np.real(qf[-1]),np.imag(qf[-1]),50,'g',marker='v',zorder=5)
    

ax[0,0].set_xlim(-1.5,1.5)
ax[0,0].set_ylim(-1.5,1.5)
ax[0,0].legend(loc=10,fontsize=9)

ax[0,0].set_ylabel(r'$v$')
ax[1,0].set_ylabel(r'$v$')
ax[1,0].set_xlabel(r'$u$')
ax[1,1].set_xlabel(r'$u$')

#plt.tight_layout()
plt.show()

    
    
    