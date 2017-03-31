## STARS

import numpy as np
import matplotlib.pyplot as plt

## parameters

m1 = 1000.
m2 = 1.
gamma = -5.
n = 1000
dt = 0.00001
e = 2

## initial conditions

r1 = np.empty((n,2))
r2 = np.empty((n,2))

r1[0,:] = np.array([0,-.1])
r2[0,:] = np.array([0,.1])

u1 = np.array([2e-1,2e-1])
u2 = np.array([-2e-1,-2e-1])

def rhs(r1,r2):
    return -gamma*(r2-r1)/((r2-r1)**2).sum()**(3./2)

# initial speed
r1[1,:] = r1[0,:] + u1*dt
r2[1,:] = r2[0,:] + u2*dt

for i in range(1,n-1):
    
    F = rhs(r1[i,:],r2[i,:])

    r1[i+1,:] = 2*r1[i,:] - r1[i-1,:] + dt**2*(m2*F - e*((r1[i,:] - r1[i-1,:])**2).sum()/dt)
    r2[i+1,:] = 2*r2[i,:] - r2[i-1,:] - dt**2*(m1*F - e*((r2[i,:] - r2[i-1,:])**2).sum()/dt)

fig1,ax1 = plt.subplots(1,1)

ax1.plot(r1[:,0],r1[:,1],color='r')
ax1.plot(r2[:,0],r2[:,1],color='g')

ax1.set_xlim(-.2,.2)
ax1.set_ylim(-.2,.2)

plt.show()