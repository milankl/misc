## THE PATTERN GENERATOR v0

import numpy as np
import time as tictoc
import matplotlib.pyplot as plt

## LOAD DATA
res = 'high' #high, low
(eofs,pcs,eigs) = np.load('python/gyres/theta_eofs_'+res+'res.npy')
print(res+' resolution data read.')

## CHOOSE MODES

modei = 500
nmode = 100
modej = modei + nmode

#mnum = np.arange(1000)
modeij = slice(modei,modej)
#modeij = (mnum < 100)+(mnum >= 900)
#nmode = np.sum(modeij)


print('Reconstruct random field from mode '+str(modei)+' to '+str(modej)+'...')

## CREATE RANDOM PCS
K = 5 # AR(K) order
exec(open('python/gyres/ar_parameters.py').read())
phi,stdnoise = ar_parameters(pcs[:,modeij],K)
print('AR'+str(K)+' parameters for '+res+' resolution computed.')

n = pcs.shape[0] # random pc length same as sample pc length
spinup = 5000 

randpc = np.zeros((nmode,n+spinup))
randpc[:,:K] = np.random.randn(nmode,K)

tic = tictoc.time()
for i in range(K-1,n-1+spinup):
    randpc[:,i+1] = (phi[:,::-1]*randpc[:,i-K+1:i+1]).sum(axis=1) \
     + stdnoise*np.random.randn(nmode)

randpc = randpc[:,spinup:] # kill spinup

## ensure mean = 0, var = 1
randpc = randpc.T - randpc.T.mean(axis=0)
randpc = (randpc / randpc.std(axis=0)).T
print('Random PCs are forced to have mean=0, var=1.')
print(str(nmode)+' PCs generated as AR'+str(K)+' process in '+str(tictoc.time()-tic)[:5]+'s.')

## MULTIPLE LINEAR REGRESSION TO PROJECT FIELD ONTO MISSING VARIANCE
c = 1.

tic = tictoc.time()
tlo =  np.load('python/gyres/temp_lowres_sfc.npy')
thi =  np.load('python/gyres/temp_highres_sfc.npy')
dvar = (thi.var(axis=0) - tlo.var(axis=0)).flatten()
print('Missing variance in '+str(tictoc.time()-tic)[:5]+'s.')

tic = tictoc.time()
esq = (eofs[modeij,...]**2).reshape((nmode,-1)).T
alpha = np.sqrt(abs(np.linalg.lstsq(esq,c*dvar)[0]))
srandpc = (randpc.T*alpha).T #scaled rand pcs
print('Multi linear regression onto missing variance in '+str(tictoc.time()-tic)[:5]+'s.')

## RECONSTRUCT FIELD
tic = tictoc.time()
randfield = np.einsum('ji,jkl',srandpc,eofs[modeij,:,:])
print('Superposition of '+str(nmode)+' EOF modes in '+str(tictoc.time()-tic)[:5]+'s.')

## SAVE

suffix = 'ar'+str(K)+'_'+str(modei)+'-'+str(modej)+'_'+res
path = 'python/gyres/patterns/randpattern_'+suffix+'.npy'

np.save(path,randfield)
print('Pattern saved in '+path)

















