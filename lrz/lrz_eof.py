## EOF TRANSFORMATION OF LORENZ SYSTEM

import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt

## define LORENZ functional F

def set_F(sig=10.,rho=28.,beta=8./3.): # LORENZ-FUNCTIONAL F
    F = np.zeros((3,10))
    F[0,1], F[0,2] = -sig, sig
    F[1,1], F[1,2], F[1,9] = rho,-1., -1.
    F[2,3], F[2,7] = -beta, 1.
    
    return F

F = set_F()

# define function H that extends xyz to its products
H = lambda x: np.array([1,x[0],x[1],x[2],x[0]**2,x[1]**2,x[2]**2,x[0]*x[1],x[1]*x[2],x[2]*x[0]])

## define right hand side of lorenz system

def rhs(xyz):
    return F.dot(H(xyz))

## parameters

dt = .001    # time spacing
N = 30000    # number of integrations

s = np.zeros((N,3))   # preallocate state vector s
# initial conditions somewhere on the attractor
s[0,:] = np.array([  2.79226311,   4.57459498,  19.59359772])    

## time integration
# using Runge Kutta 4 scheme

for i in range(N-1):
    k1 = rhs(s[i,:])
    k2 = rhs(s[i,:] + .5*dt*k1)
    k3 = rhs(s[i,:] + .5*dt*k2)
    k4 = rhs(s[i,:] + dt*k3)
    s[i+1,:] = s[i,:] + dt/6. * (k1 + 2*k2 + 2*k3 + k4)

## EOF transformation
solver = Eof(s)

##
# scaling either 0,1 or 2
e = solver.eofs(eofscaling=1)

eF = e.T.dot(F)
sm = s.mean(axis=0)

def rhs_eof(abc):
    #abc is the state vector in EOF coordinates
    return eF.dot(H(abc.dot(e.T) + sm))

## RUNGE KUTTA in EOF coordinates
dt = .5
N = 10000

a = np.empty((N,3))   # preallocate state vector a
# initial conditions somewhere on the attractor
a[0,:] = np.array([ 8.04669874,  7.76236837, -3.39516727])

for i in range(N-1):
    k1 = rhs_eof(a[i,:])
    k2 = rhs_eof(a[i,:] + .5*dt*k1)
    k3 = rhs_eof(a[i,:] + .5*dt*k2)
    k4 = rhs_eof(a[i,:] + dt*k3)
    a[i+1,:] = a[i,:] + dt/6. * (k1 + 2*k2 + 2*k3 + k4)

##

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(a[:,0],a[:,1],a[:,2])

plt.show()



"""
m = solver.eigenvalues()
p = solver.pcs(pcscaling=0)

r = solver.reconstructedField(neofs=3) + s.mean(axis=0)
q = p.dot(e) + s.mean(axis=0)

fig,axs = plt.subplots(1,3,sharex=True)

for j in range(3):
    axs[j].plot(s[:,j],lw=2)
    axs[j].plot(r[:,j],'w--')
    axs[j].plot(q[:,j],'k-.')

plt.show()
"""