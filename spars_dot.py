import numpy as np
from scipy import sparse

n = 25000

A = sparse.diags(np.random.randint(0,10,n),0,shape=(n,n)).tocsr() \
    + sparse.diags(np.random.randint(0,10,n-1),1,shape=(n,n)).tocsr()

x = np.random.randint(0,10,n)

Ad = A.todense()

# manually

def sdot(A,x):
    b = np.empty_like(x)
    for i in range(n):
        u,v = A.indptr[i:i+2]
        a = A.data[u:v]
        c = x[A.indices[u:v]]    
        b[i] = sum(a*c)
    
    return b