## EOF TRANSFORMATION OF LORENZ SYSTEM

import numpy as np
from eofs.standard import Eof
import matplotlib.pyplot as plt

## define right hand side of lorenz system

def rhs(xyz,sig=10.,rho=28.,beta=8./3.):
    xdot = sig*(xyz[1] - xyz[0])
    ydot = xyz[0]*(rho - xyz[2]) - xyz[1]
    zdot = xyz[0]*xyz[1] - beta*xyz[2]
    return np.array([xdot,ydot,zdot])

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
# scaling either (0,0), (1,2) or (2,1)
e = solver.eofs(eofscaling=0)

def set_F(e,sm,sig=10.,rho=28.,beta=8./3.): # LORENZ-EOF-FUNCTIONAL F
    # G1,G2,G3 are the columns of the functional xyz_dot = G(a1,a2,a3)
    # e is the orthogonal EOF matrix
    # sm is the time mean of the state vector s
    G1 = [sig*(sm[1]-sm[0]),sig*(e[0,1]-e[0,0]),sig*(e[1,1]-e[1,0]),\
        sig*(e[2,1] - e[2,0]),0,0,0,0,0,0]
    G2 = [rho*(sm[1]-sm[0])-sm[0]*sm[2],\
        rho*(e[0,0] - e[0,1]) - (sm[0]*e[0,2] + sm[2]*e[0,0]),\
        rho*(e[1,0] - e[1,1]) - (sm[0]*e[1,2] + sm[2]*e[1,0]),\
        rho*(e[2,0] - e[2,1]) - (sm[0]*e[2,2] + sm[2]*e[2,0]),\
        -e[0,0]*e[0,2],-e[1,0]*e[1,2],-e[2,0]*e[2,2],\
        e[0,0]*e[1,2] + e[1,0]*e[0,2],\
        e[1,0]*e[2,2] + e[2,0]*e[1,2],\
        e[2,0]*e[0,2] + e[0,0]*e[2,2]]
    G3 = [sm[0]*sm[1]-beta*sm[2],\
        sm[0]*e[0,1] + sm[1]*e[0,0] - beta*e[0,2],\
        sm[0]*e[1,1] + sm[1]*e[1,0] - beta*e[1,2],\
        sm[0]*e[2,1] + sm[1]*e[2,0] - beta*e[2,2],\
        e[0,0]*e[0,1],e[1,0]*e[1,1],e[2,0]*e[2,1],\
        e[0,0]*e[1,1] + e[1,0]*e[0,1],\
        e[1,0]*e[2,1] + e[2,0]*e[1,1],\
        e[2,0]*e[0,1] + e[0,0]*e[2,1]]
    
    # return F, as a_i,dot = F(a_i) = e_inv * G(a_i)
    return e.T.dot(np.array((G1,G2,G3)))

F = set_F(e,s.mean(axis=0))

def rhs_eof(abc):
    #abc is the state vector in EOF coordinates
    #F is taken from outside the function
    ABC = np.array([1,abc[0],abc[1],abc[2],\
        abc[0]**2,abc[1]**2,abc[2]**2,\
        abc[0]*abc[1],abc[1]*abc[2],abc[2]*abc[0]])
    return F.dot(ABC)

## RUNGE KUTTA in EOF coordinates
a = np.empty((N,3))   # preallocate state vector a
# initial conditions somewhere on the attractor
a[0,:] = np.array([.01, .01,.01])

dt = 1e-4
N = 1000

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