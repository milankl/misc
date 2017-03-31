import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

H = 1.
nz = 100
dz = H/nz

#initial conditions
a = np.zeros(nz)
#a[30:60] = 1

# laplace operator
L = (sparse.diags(-2*np.ones(nz),0,shape=(nz,nz)) +\
    sparse.diags(np.ones(nz),1,shape=(nz,nz)) +\
    sparse.diags(np.ones(nz),-1,shape=(nz,nz))).tocsr() / dz**2

#L[0,0] = -1/dz**2
#L[-1,-1] = -1/dz**2
L = L.tocsr()

nu = 1
dt = dz**2/nu*0.5

M = sparse.eye(nz) - nu*dt*L

s = np.zeros(nz)
s[0] = 1/dz**2*nu

for i in range(100):
    a = spsolve(M,a+dt*s)

print(a.mean())    
plt.plot(np.hstack((1,a,0)))
plt.show()