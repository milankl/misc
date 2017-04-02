## AIC/BIC criterion

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc
from scipy.stats import rankdata
exec(open('python/ecco2/local_functions.py').read())
exec(open('python/ecco2/colormaps.py').read())

#(eofs,pcs,eigs) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs,pcs,eigs) = np.load('python/gyres/theta_eofs_highres.npy')

modmax = 100 #loop until this mode
K = 20 #MAXIMUM ORDER
n = pcs.shape[0] #length of timeseries

aic = np.zeros((K+1,modmax))
bic = np.zeros((K+1,modmax))

for m in range(modmax): #loop over modes
    r = acf(pcs[:,m],K+1)[1:] # sample autocorrelation

    a = np.zeros((K,K))
    a[0,0] = r[0]

    for k in range(2,K+1):
        A = np.eye(k)
        for i in range(1,k):
            A = A + np.diagflat([r[i-1]]*(k-i),i) + np.diagflat([r[i-1]]*(k-i),-i)
        a[k-1,:k] = np.linalg.solve(A,r[:k])

    snoise = np.hstack((1,np.sqrt(abs(1 - np.dot(a,r)))))
    ordvec = np.arange(K+1)
    aic[:,m] = n*np.log(n/(n-ordvec-1)*snoise) + 2*(ordvec + 1)
    bic[:,m] = n*np.log(n/(n-ordvec-1)*snoise) + (ordvec + 1)*np.log(n)

## VISUALIZE THE RESULT

aicmin = np.nanargmin(aic,axis=0)
bicmin = np.nanargmin(bic,axis=0)

aicrank = np.array([rankdata(aic[:,m]) for m in range(modmax)])
bicrank = np.array([rankdata(bic[:,m]) for m in range(modmax)])

## PLOTTING

plt.figure(1)

x = np.arange(modmax+1) - .5
y = np.arange(K+2) - .5

xx,yy = np.meshgrid(x,y)

plt.pcolor(xx,yy,bicrank.T,cmap=viridis_r)
plt.plot(aicmin,'ko',label='aic')
plt.plot(bicmin,'wo',label='bic')
plt.xlim(-.5,modmax-.5)
plt.ylim(-.5,K+.5)
plt.xlabel('mode #')
plt.ylabel('AR order')
plt.title('BIC rank')
plt.legend(loc=1)
plt.colorbar()

plt.show()
