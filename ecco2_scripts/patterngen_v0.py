## THE PATTERN GENERATOR v0

import numpy as np
import time as tictoc
import matplotlib.pyplot as plt

## LOAD DATA

time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(eofs,pcs,eigs) = np.load('python/ecco2/theta_eofs_1000.npy')

## CHOOSE MODES

modei = 50
modej = modei + 100
nmode = modej - modei

print('Reconstruct random field from mode '+str(modei)+' to '+str(modej)+'...')

## CREATE RANDOM PCS

phi = np.load('python/ecco2/ar4_parameters.npy')
stdnoise = phi[modei:modej,-1]
phi = phi[modei:modej,:-1]

K = phi.shape[1] # AR(K) order
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

## RECONSTRUCT FIELD
tic = tictoc.time()
randfield = np.einsum('ji,jkl',randpc,eofs[modei:modej,:,:])
print('Superposition of '+str(nmode)+' EOF modes in '+str(tictoc.time()-tic)[:5]+'s.')

## MASKING
tic = tictoc.time()
randfield = np.ma.masked_array(randfield,mask=np.array([eccomask_NA]*n))
print('Masking in '+str(tictoc.time()-tic)[:5]+'s.')

## SAVE

path = 'python/ecco2/patterns/randpattern_4.npy'

np.save(path,(randfield,[modei,modej]))
print('Pattern saved in '+path)

















