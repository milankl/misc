## SAVE AUTOREGRESSIVE PROCESS PARAMETERS

def ar_parameters(pcs,K):

    import numpy as np
    exec(open('python/ecco2/local_functions.py').read())

    # the AR parameters and standard deviation of white noise
    phi = np.zeros((pcs.shape[1],K))
    stdnoise = np.zeros(pcs.shape[1]) 

    for m in range(pcs.shape[1]): # loop over modes
        
        r = acf(pcs[:,m],K+1)[1:] # sample autocorrelation function

        # set up Yule-Walker-Matrix A
        A = np.eye(K)
        for i in range(1,K):
            A = A + np.diagflat([r[i-1]]*(K-i),i)
        A = A + A.T - np.eye(K)

        # solve Yule-Walker-equation system
        phi[m,:] = np.linalg.solve(A,r)

        # standard deviation for white noise
        stdnoise[m] = np.sqrt(abs(1 - np.dot(phi[m,:K],r)))
    
    return phi,stdnoise
    
