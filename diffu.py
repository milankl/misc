import numpy as np
import matplotlib.pyplot as plt

N = 300

x = np.linspace(-2,2,N)
dx = x[1]-x[0]

##
u = 2*(x > 0) - 1

#h = np.exp(-x**2*4)
#u = np.gradient(h)
##

L = (np.diag(-2*np.ones(N),0) + np.diag(np.ones(N-1),1) + np.diag(np.ones(N-1),-1)) / dx**2
L2 = L.dot(L)

nu = 2
dt = 0.00001

uh = u.copy()
ub = u.copy()
ue = u.copy()

for i in range(1000):
    uh = uh + dt*nu*L.dot(uh)
    ub = ub - dt*nu*dx**2*L2.dot(ub) 
    ue = ue - dt*nu*dx**2*L2.dot(ue) - 0.1*dt*nu*L.dot(ue)

##
fig,ax = plt.subplots(1,1)

#ax.plot(x,u,'-',drawstyle='steps-mid',label=r'$u_0$')
#ax.plot(x,uh,'-',drawstyle='steps-mid',label=r'$\partial_tu = \nu_A\partial_x^2u$')
#ax.plot(x,ub,'-',drawstyle='steps-mid',label=r'$\partial_tu = -\nu_B\partial_x^4u$')
#ax.plot(x,ue,'-',drawstyle='steps-mid',label=r'$\partial_tu = -\nu_B\partial_x^4u - \nu_{back}\partial_x^2u$')

ax.plot(x,u,'-',label=r'$u_0$')
ax.plot(x,uh,'-',label=r'$\partial_tu = \nu_A\partial_x^2u$')
ax.plot(x,ub,'-',label=r'$\partial_tu = -\nu_B\partial_x^4u$')
ax.plot(x,ue,'-',label=r'$\partial_tu = -\nu_B\partial_x^4u + \nu_{back}\partial_x^2u$')

ax.plot(x,np.zeros_like(x),'--',color='grey')

ax.set_xlabel('x')
ax.set_ylabel('u')

ax.set_xlim(-1,1)
ax.legend(loc=2)


plt.show()