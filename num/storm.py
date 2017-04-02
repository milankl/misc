## MOVING STORM
import numpy as np
import matplotlib.pyplot as plt
import pylan as pn

##grid
#in meters
dx = 25e3
dy = 25e3

x = np.arange(-1e6,1e6+dx,dx)
y = np.arange(-1e6,1e6+dy,dy)

xx,yy = np.meshgrid(x,y)

#parameters
f = 1e-4
p0 = 200 #in Pascal
L = 1e5
y0 = -5e5
vs = 5 * 1e3/3600. #meter per second

#time
T = 1e6 / vs

## pressure p, forcing X

def p(t):
    r2 = xx**2 + (yy - y0 - vs*t)**2
    return -p0*np.exp(-r2/L**2)

def X(p):
    dpdx,dpdy = np.gradient(p.T,dx,dy)
    return -(dpdx.T + 1j*dpdy.T)

## euler forward
# initial conditions
q0 = X(p(0))/1j/f

def q(dt):
    n = int(T/dt)
    t = np.arange(0,T,dt)
    q = np.zeros((n+1,)+q0.shape,dtype=complex)
    q[0,...] = q0
    #Euler forward for first step
    q[1,...] = (1 - 1j*f*dt)*q[0,...] + dt*X(p(t[0]))
    
    #Leap-frog for the rest
    for i in range(1,n):
        q[i+1,...] = q[i-1,...] + 2*dt*(X(p(t[i])) - 1j*f*q[i,...])
    
    return q

#solution for q
q1 = q(1000)

#indices
n1 = int(q1.shape[0]/4-.25)*np.arange(5)

#pressure
pp = np.array([p(ti) for ti in T/4*np.arange(5)])

## plotting

fig,axs = plt.subplots(1,5,sharex=True,sharey=True)
qks = [0,0,0,0,0]

for i,i1 in zip(range(5),n1):
    
    im1 = axs[i].contourf(xx/1e5,yy/1e5,pp[i]/100)
    qks[i] = axs[i].quiver(xx/1e5,yy/1e5,np.real(q1[i1,...]),np.imag(q1[i1,...]),scale=150,width=0.012)
    axs[i].set_title(r'$t$ = %ih' % (T/4*np.arange(5)/3600)[i])
    

axs[0].set_xlim(-1.8*L/1e5,1.8*L/1e5)
axs[0].set_ylim((y0-2*L)/1e5,-(y0-2*L)/1e5)
for i in range(5):
    axs[i].set_xticks([-1,0,1])
    
axs[2].set_xlabel(r'$x$ / 100km')
axs[0].set_ylabel(r'$y$ / 100km')
cb = plt.colorbar(im1,ax=(axs[0],axs[1],axs[2],axs[3],axs[4]))
cb.set_label('pressure [hPa]')

axs[0].add_patch(plt.Rectangle([-1.3,4.5],2.6,1.4,color='white'))
axs[0].add_patch(plt.Rectangle([-1.3,4.5],2.6,1.4,fill=False,color='k'))
axs[0].quiverkey(qks[0], 0.5, 0.85, 30, '30 m/s')

plt.show()

    