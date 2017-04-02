## SAVE AUTOREGRESSIVE PROCESS PARAMETERS

import numpy as np
exec(open('python/ecco2/local_functions.py').read())

## LOAD PRINCIPAL COMPONENTS
(eofs,pcs,eigs) = np.load('python/ecco2/theta_eofs_1000.npy')

K = 4 #AR(k) order

# the AR parameters and standard deviation of white noise
phi = np.zeros((pcs.shape[1],K+1)) 

for m in range(pcs.shape[1]): # loop over modes
        
    r = acf(pcs[:,m],K+1)[1:] # sample autocorrelation function

    # set up Yule-Walker-Matrix A
    A = np.eye(K)
    for i in range(1,K):
        A = A + np.diagflat([r[i-1]]*(K-i),i)
    A = A + A.T - np.eye(K)
    
    # solve Yule-Walker-equation system
    phi[m,:K] = np.linalg.solve(A,r)

    # standard deviation for white noise
    phi[m,-1] = np.sqrt(1 - np.dot(phi[m,:K],r))

np.save('python/ecco2/ar'+str(K)+'_parameters.npy',phi)
print('AR'+str(K)+' parameters written to file.')
