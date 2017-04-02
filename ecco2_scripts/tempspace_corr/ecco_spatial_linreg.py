## ECCO LINEAR REGRESSION

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc
exec(open('python/ecco2/colormap.py').read())

## LOAD DATA

#theta = np.load('python/ecco2/theta_sfc_NA_ano_detr.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

#theta = np.ma.masked_array(theta,mask=np.array([eccomask_NA]*len(time)))
#theta = theta - theta.mean(axis=0)

## PIC ONE POINT

lat1 = np.argmin(abs(latna - 48.3))
lon1 = np.argmin(abs(lonna - (323.15)))
theta1 = theta[:,lat1,lon1]

## DO SPATIAL CORRELATION

def mlinreg(x,y):
    """linear regression slope between vector x and 3d-array y. faster than ployfit. x,y are assumed to have zero mean in time dimension (first dim)."""
    yr = np.ma.reshape(y,(y.shape[0],np.prod(y.shape[1:])))
    return (((np.ma.dot(x,yr).reshape(y.shape[1:])) / (x.shape[0] - 1)) / x.var())

tic = tictoc.time()
scorr = mlinreg(theta1,theta)
print(tictoc.time() - tic)

## PLOTTING

plt.pcolormesh(lonna,latna,scorr,cmap=viridis)
plt.colorbar()
#plt.clim(-1,1)
plt.contour(lonna,latna,scorr,np.ones((3,))*.3,colors='k')
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.show()

