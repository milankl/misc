## WAVE EQ
# solving dudt = -gdhdx & dhdt = -dx(Hu)
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

##

nx = 400
nt = 800
Lx = 1.
sigma=0.08
cfl = 1.5
g = 20.

#depth
H = sparse.csc_matrix(np.diag(np.ones(nx)))

#domain
dx = Lx/nx
u = np.sqrt(g*H.max())
dt = cfl * dx / u

#grid
x = np.arange(nx)*dx

#indices for producing gradient matrices
p = np.arange(nx)
pm = (p-1)%nx
pp = (p+1)%nx

##2nd order gradient
i = np.concatenate((p,p))
j = np.concatenate((p,pp))
s = np.concatenate((p*0-1,p*0+1))/(2*dx)

# von neumann
s[-1] = 0
s[nx-1] = 0

# dirichlet
#s[-1] = 0
#s[0] = 0

G2x = sparse.coo_matrix((s,(i,j)),shape=(nx,nx)).tocsc()

##initial state
#s = np.zeros(nx)
s = np.exp(-(x-.4)**2/(2*sigma**2))#*(x>.3)*(x<.5)*1

M = -u*G2x
M = -(M.T.dot(H)).dot(M)
q = (1-cfl)/(1+cfl)

##
fig = plt.figure()
ax = fig.add_subplot(111)
sout = np.zeros((nx,nt))

for kt in range(nt):
    
    sout[:,kt] = s
    ds = dt**2*M.dot(s)
    
    if kt == 0:
        s = s + ds
    else:
        s = 2*s - sout[:,kt-1] + ds
        s[-1] = sout[-2,kt-1] - q*s[-2] + q*sout[-1,kt-1]

    if kt == 0:
        l1, = plt.plot(x[1:],s[1:],lw=2)
        plt.ylim(-1.5,1.5)
    else:
        plt.pause(0.00001)
        l1.set_data(x[1:],s[1:])

plt.pause(1)
plt.close('all')