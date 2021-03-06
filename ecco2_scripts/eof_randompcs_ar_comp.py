## PRODUCE RANDOM PRINCIPAL COMPONENTS

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/local_functions.py').read())

(eofL,pctau) = np.load('python/ecco2/eof_ltscales.npy')
(thetaeofs,thetapcs,eigs) = np.load('python/ecco2/theta_eofs_1000.npy')

mode = 12
K = 1 #AR(k) order

r = acf(thetapcs[:,mode],K+1)[1:]

## via system of equations

a = np.zeros((K,K))
a[0,0] = r[0]

for k in range(2,K+1):
    A = np.eye(k)
    for i in range(1,k):
        A = A + np.diagflat([r[i-1]]*(k-i),i) + np.diagflat([r[i-1]]*(k-i),-i)
    a[k-1,:k] = np.linalg.solve(A,r[:k])

snoise = np.sqrt(1 - np.dot(a,r))

n = 100
l = 60
spinup = 1000
acfs = np.zeros((n,l))

for m in range(n):

    randpc = np.zeros(thetapcs[:,0].shape[0]+spinup)
    randpc[:K] = np.random.randn(K)

    for i in range(K-1,thetapcs[:,0].shape[0]-1+spinup):
        randpc[i+1] = np.dot(a[K-1,::-1],randpc[i-K+1:i+1]) \
         + snoise[K-1]*np.random.randn(1)
    
    randpc = randpc[spinup:]

    acfs[m,:] = acf(randpc,l)

## --------------------------------

fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(10,10))

ax1.errorbar(np.arange(l)*3,acfs.mean(axis=0),yerr=2*acfs.std(axis=0))
ax1.plot(np.arange(l)*3,acf(thetapcs[:,mode],l),'r')
ax1.plot(np.arange(l)*3,np.ones(l)/np.exp(1),'grey')
ax1.set_ylim(-.4,1)
ax1.set_title('Pot. Temp: PC of mode 12, AR(1)-Process')

## ------------------------------------

mode = 12
K = 4 #AR(k) order

r = acf(thetapcs[:,mode],K+1)[1:]

## via system of equations

a = np.zeros((K,K))
a[0,0] = r[0]

for k in range(2,K+1):
    A = np.eye(k)
    for i in range(1,k):
        A = A + np.diagflat([r[i-1]]*(k-i),i) + np.diagflat([r[i-1]]*(k-i),-i)
    a[k-1,:k] = np.linalg.solve(A,r[:k])

snoise = np.sqrt(1 - np.dot(a,r))

n = 100
l = 60
spinup = 1000
acfs = np.zeros((n,l))

for m in range(n):

    randpc = np.zeros(thetapcs[:,0].shape[0]+spinup)
    randpc[:K] = np.random.randn(K)

    for i in range(K-1,thetapcs[:,0].shape[0]-1+spinup):
        randpc[i+1] = np.dot(a[K-1,::-1],randpc[i-K+1:i+1]) \
         + snoise[K-1]*np.random.randn(1)
    
    randpc = randpc[spinup:]

    acfs[m,:] = acf(randpc,l)

ax2.errorbar(np.arange(l)*3,acfs.mean(axis=0),yerr=2*acfs.std(axis=0))
ax2.plot(np.arange(l)*3,acf(thetapcs[:,mode],l),'r')
ax2.plot(np.arange(l)*3,np.ones(l)/np.exp(1),'grey')
ax2.set_ylim(-.4,1)
ax2.set_title('..., AR(4)-Process')

## ------------------------------------

mode = 12
K = 12 #AR(k) order

r = acf(thetapcs[:,mode],K+1)[1:]

## via system of equations

a = np.zeros((K,K))
a[0,0] = r[0]

for k in range(2,K+1):
    A = np.eye(k)
    for i in range(1,k):
        A = A + np.diagflat([r[i-1]]*(k-i),i) + np.diagflat([r[i-1]]*(k-i),-i)
    a[k-1,:k] = np.linalg.solve(A,r[:k])

snoise = np.sqrt(1 - np.dot(a,r))

n = 100
l = 60
spinup = 1000
acfs = np.zeros((n,l))

for m in range(n):

    randpc = np.zeros(thetapcs[:,0].shape[0]+spinup)
    randpc[:K] = np.random.randn(K)

    for i in range(K-1,thetapcs[:,0].shape[0]-1+spinup):
        randpc[i+1] = np.dot(a[K-1,::-1],randpc[i-K+1:i+1]) \
         + snoise[K-1]*np.random.randn(1)
    
    randpc = randpc[spinup:]

    acfs[m,:] = acf(randpc,l)

ax3.errorbar(np.arange(l)*3,acfs.mean(axis=0),yerr=2*acfs.std(axis=0),label='AR mean +- 2std')
ax3.plot(np.arange(l)*3,acf(thetapcs[:,mode],l),'r',label='observed acf')
ax3.plot(np.arange(l)*3,np.ones(l)/np.exp(1),'grey',label='threshold')
ax3.set_ylim(-.4,1)
ax3.legend(loc=1)
ax3.set_title('..., AR(12)-Process')
ax3.set_xlabel('lag [days]')

plt.show()

