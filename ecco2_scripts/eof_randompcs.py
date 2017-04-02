## PRODUCE RANDOM PRINCIPAL COMPONENTS
## with AR1

import numpy as np
import matplotlib.pyplot as plt

def acf(x,l):
	return np.array([1]+[np.corrcoef(x[:-i],x[i:])[0,1] for i in range(1,l)])

(eofL,pctau) = np.load('python/ecco2/eof_ltscales.npy')
(thetaeofs,thetapcs,eigs) = np.load('python/ecco2/theta_eofs_1000.npy')

mode = 97
(r1,r2) = acf(thetapcs[:,mode],3)[1:]

alpha1 = r1*(1-r2)/(1-r1**2)
alpha2 = (r2 - r1**2) / (1-r1**2)

## via system of equations

A = np.array([[1, r1],[r1, 1]])
(alpha1la,alpha2la) = np.linalg.solve(A,np.array([r1,r2]))

n = 50
l = 100
acfs = np.zeros((n,l))

for m in range(n):

    randpc = np.zeros(thetapcs[:,0].shape[0]+1000)
    randpc[:2] = np.random.randn(2)

    for i in range(1,thetapcs[:,0].shape[0]-1+1000):
        randpc[i+1] = alpha1*randpc[i] + alpha2*randpc[i-1] + np.sqrt((1-alpha2**2)*(1-r1**2))*np.random.randn(1)
        
    randpc = randpc[-thetapcs[:,0].shape[0]:]

    acfs[m,:] = acf(randpc,l)


fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(10,10))

ax1.plot(thetapcs[:,mode],'r')
ax1.set_title('std = '+str(np.std(thetapcs[:,mode]))[:5])

ax2.plot(randpc)
ax2.set_title('std = '+str(np.std(randpc))[:5])

ax3.errorbar(np.arange(l),acfs.mean(axis=0),yerr=2*acfs.std(axis=0))
ax3.plot(np.arange(l),acf(thetapcs[:,mode],l),'r')
ax3.plot(np.arange(l),np.ones(l)/np.exp(1),'grey')

plt.show()

