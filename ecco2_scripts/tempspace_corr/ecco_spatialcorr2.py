## ECCO SPATIAL CORRELATION

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc
import scipy.stats as stats


## LOAD DATA

theta = np.load('python/ecco2_dT/dtheta_sfc_NA_ano_v2.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

## REMOVE MEAN
tic = tictoc.time()
theta = theta - theta.mean(axis=0)
theta = theta / theta.std(axis=0)
print('Normalise, standardise in '+str(tictoc.time()-tic)[:5]+'s.')

## SPATIAL CORRELATION FUNCTION
def mcorr2(x,y):
    """correlation between vector x and 2d-array y. faster than np.corrcoef. x,y are assumed to have zero mean in time dimension (first dim) and variance = 1."""
    yr = np.ma.reshape(y,(y.shape[0],np.prod(y.shape[1:])))
    return (np.ma.dot(x,yr).reshape(y.shape[1:]) / (x.shape[0] - 1))

## PIC ONE POINT

lat1 = np.argmin(abs(latna - 30.))
lon1 = np.argmin(abs(lonna - (285.)))
theta1 = theta[:,lat1,lon1]

## DO SPATIAL CORRELATION

def mcorr(x,y):
    """correlation between vector x and 3d-array y. faster than np.corrcoef. x,y are assumed to have zero mean in time dimension (first dim)."""
    yr = np.ma.reshape(y,(y.shape[0],np.prod(y.shape[1:])))
    return (((np.ma.dot(x,yr).reshape(y.shape[1:])) / (x.shape[0] - 1) / y.std(axis=0)) / x.std())
    
def rcorr(x,y):
    """rank correlation"""
    tic = tictoc.time()
    x = stats.rankdata(x)
    print(tictoc.time() - tic)
    tic = tictoc.time()
    y = np.rollaxis(np.apply_along_axis(stats.rankdata,0,y),0,len(y.shape))
    print(tictoc.time() - tic)
    tic = tictoc.time()
    z = (1 - 6*((y-x)**2).sum(axis=(len(y.shape)-1))/x.shape[0]/(x.shape[0]**2-1))
    print(tictoc.time() - tic)
    return z

tic = tictoc.time()
scorr = mcorr2(theta1,theta)
print(tictoc.time() - tic)

tic = tictoc.time()
#scorr2 = rcorr(theta1,theta)
print(tictoc.time() - tic)


## PLOTTING

plt.figure(1)
plt.pcolormesh(lonna,latna,scorr)
plt.colorbar()
plt.clim(-1,1)
plt.contour(lonna,latna,scorr,np.ones((3,))/np.exp(1),colors='k')
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))

#plt.figure(2)
#plt.pcolormesh(lonna,latna,scorr2)
#plt.colorbar()
#plt.clim(-1,1)
#plt.contour(lonna,latna,scorr,np.ones((3,))/np.exp(1),colors='k')
#plt.xlim((lonna[0],lonna[-1]))
#plt.ylim((latna[0],latna[-1]))
plt.show()

