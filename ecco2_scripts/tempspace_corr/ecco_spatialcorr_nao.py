## ECCO SPATIAL CORRELATION

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc
import scipy.stats as stats
exec(open('python/ecco2/colormap.py').read())


## LOAD DATA

#theta = np.load('python/ecco2/theta_sfc_NA_ano_detr.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

## LOAD NAO

dat = np.loadtxt('python/nao/norm.daily.nao.index.b500101.current.ascii')
naoyear = dat[:,0]
naomon = dat[:,1]
naoday = dat[:,2]
nao = dat[:,3]

nao[nao == -999.9] = None
ibegin = np.where((naoyear == 1992) * (naomon == 1) * (naoday == 1))[0][0]
iend = ibegin + 3*len(time)

nao = np.nanmean(nao[ibegin:iend].reshape(3,-1,order='F'),axis=0)
nao = -(nao - nao.mean())
nao = nao / nao.std()
## REMOVE MEAN

theta = theta - theta.mean(axis=0)

## PIC ONE POINT

lat1 = np.argmin(abs(latna - 15.))
lon1 = np.argmin(abs(lonna - (325.)))
theta1 = theta[:,lat1,lon1]

## DO SPATIAL CORRELATION

def mcorr(x,y):
    #correlation between vector x and 3d-array y. faster than np.corrcoef. x,y are assumed to have zero mean in time dimension (first dim).
    yr = np.ma.reshape(y,(y.shape[0],np.prod(y.shape[1:])))
    return (((np.ma.dot(x,yr).reshape(y.shape[1:])) / (x.shape[0] - 1) / y.std(axis=0)) / x.std())

def linreg(x,y):
    #linear regression between vector x and 3d-array y. x,y are assumed to have zero mean in time dimension (first dim).
    yr = np.ma.reshape(y,(y.shape[0],np.prod(y.shape[1:])))
    return (((np.ma.dot(x,yr).reshape(y.shape[1:])) / (x.shape[0] - 1)) / x.var())   


tic = tictoc.time()
scorr = mcorr(theta1,theta)
print(tictoc.time() - tic)

tic = tictoc.time()
regr = linreg(theta1,theta)
print(tictoc.time() - tic)

tic = tictoc.time()
naoscorr = mcorr(nao,theta)
print(tictoc.time() - tic)

tic = tictoc.time()
naoregr = linreg(nao,theta)
print(tictoc.time() - tic)

scorr.mask = eccomask_NA
regr.mask = eccomask_NA
naoscorr.mask = eccomask_NA
naoregr.mask = eccomask_NA


## PLOTTING

xx,yy = np.meshgrid(lonna,latna)

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(16,10))

im1 = ax1.pcolormesh(xx,yy,scorr)
ax1.contour(xx,yy,scorr,np.ones((3,))/np.exp(1),colors='k')
ax1.set_xlim((lonna[0],lonna[-1]))
ax1.set_ylim((latna[0],latna[-1]))
ax1.set_title('Pot. Temp: Correlation with 15N, 325E')
plt.colorbar(im1,ax=ax1)
im1.set_clim(-.5,1)

im2 = ax2.pcolormesh(xx,yy,regr)
ax2.set_xlim((lonna[0],lonna[-1]))
ax2.set_ylim((latna[0],latna[-1]))
ax2.set_title('Pot. Temp: Lin regr with 15N, 325E')
plt.colorbar(im2,ax=ax2)
im2.set_clim(-.5,1)

im3 = ax3.pcolormesh(xx,yy,naoscorr)
ax3.set_xlim((lonna[0],lonna[-1]))
ax3.set_ylim((latna[0],latna[-1]))
ax3.set_title('Pot. Temp: Correlation with NAO')
plt.colorbar(im3,ax=ax3)
im3.set_clim(-.2,.2)

im4 = ax4.pcolormesh(xx,yy,naoregr)
ax4.set_xlim((lonna[0],lonna[-1]))
ax4.set_ylim((latna[0],latna[-1]))
ax4.set_title('Pot. Temp: Lin regr with NAO')
plt.colorbar(im4,ax=ax4)
im4.set_clim(-.2,.2)

plt.show()

