## ECCO SPATIAL CORRELATION

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc

## LOAD DATA

theta = np.load('python/ecco2/theta_sfc_NA_ano_detr.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

## PIC ONE POINT

lat1 = np.argmin(abs(latna - 32.))
lon1 = np.argmin(abs(lonna - (360.-78.)))
theta1 = theta[:,lat1,lon1]

## DO SPATIAL CORRELATION

tic = tictoc.time()

scorr = np.zeros(theta[0,:,:].shape)

#APPLY MASK
scorr = np.ma.masked_array(scorr,mask=eccomask_NA)

for ilon in range(lonna.shape[0]):
    if ilon%30 == 0:
        print('processing '+str(round(ilon*1.0/(lonna.shape[0])*100))+'%')
    
    for ilat in range(latna.shape[0]):
        if ~eccomask_NA[ilat,ilon]:
            scorr[ilat,ilon] = np.corrcoef(theta1,theta[:,ilat,ilon])[0,1]

print(tictoc.time() - tic)

## PLOTTING

plt.pcolormesh(lonna,latna,scorr)
plt.colorbar()
plt.contour(lonna,latna,scorr,np.ones((3,))*.3,colors='k')
plt.clim((-1,1))
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.show()
