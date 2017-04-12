# cellular automaton model

import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse

## initial
n = 100
m = 1000
s = np.zeros((n**2),dtype=int)
s[5550] = 1
s[5551] = 1
s[5552] = 1
s[5547] = 1
s[5546] = 1
s[5549-n] = 1
s[5547-2*n] = 1


## operator matrix

D = np.zeros((n**2,n**2),dtype=int)
i = np.arange(n**2)

D[i,(i+1)%n**2] = 1
D[i,(i-1)] = 1
D[i,(i+n)%n**2] = 1
D[i,(i-n)] = 1

D[i,(i+n+1)%n**2] = 1
D[i,(i+n-1)%n**2] = 1

D[i,(i-n-1)] = 1
D[i,(i-n+1)] = 1

D = sparse.csr_matrix(D)

##

fig, ax = plt.subplots()

Sxy = s.reshape((n,n))
h = ax.spy(Sxy)
plt.pause(2)

for j in range(m):
    ds = D.dot(s)
    dk3 = abs(np.sign(ds-3))
    dk2 = abs(np.sign(ds-2))
    
    s = np.sign(s*(2-dk3-dk2) + (1-dk3))
    Sxy = s.reshape((n,n))
    h.set_data(Sxy)
    plt.pause(0.1)


plt.pause(5)
plt.close('all')

##

#x 0 1 2 3 4 5 6 7 8
#0 0 0 0 1 0 0 0 0 0
#1 0 0 1 1 0 0 0 0 0

