## VARIANCE OF HIGH VS LOW

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormap.py').read())
import scipy.stats as stats

## load data

thi = np.load('python/gyres/temp_highres_sfc.npy')
tlo = np.load('python/gyres/temp_lowres_sfc.npy')
randfield = np.load('python/gyres/patterns/randpattern_ar5_0-100+900-1000_high.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

## variances
varhi = thi.var(axis=0)
varlo = tlo.var(axis=0)
varrand = randfield.var(axis=0)

"""
tlodt,tlody,tlodx = np.gradient(tlo)

x1 = abs(tlodt).mean(axis=0).flatten()
x2 = abs(tlodx).mean(axis=0).flatten()
x3 = tlodt.var(axis=0).flatten()
x4 = varlo.flatten()
x5 = stats.skew(tlo,axis=0).flatten()
x6 = stats.kurtosis(tlo,axis=0).flatten()

X = np.vstack((x1,x2,x3,x4,x5,x6))
X = X.T / X.std(axis=1)

b = np.linalg.lstsq(X,varhi.flatten())[0]

varest = np.dot(b,X.T).reshape((128,128))
"""

##
plt.figure(1)
plt.pcolormesh(lon,lat,varhi)
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
plt.colorbar()
plt.clim(0,8)
plt.title('Var(T_high)')

plt.figure(2)
plt.pcolormesh(lon,lat,varlo)
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
plt.colorbar()
plt.clim(0,8)
plt.title('Var(T_low)')


plt.figure(3)
plt.pcolormesh(lon,lat,varhi-varlo)
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
plt.colorbar()
plt.clim(0,8)
plt.title('Var(T_high) - Var(T_low)')

plt.figure(4)
plt.pcolormesh(lon,lat,varrand)
plt.xlim(lon.min(),lon.max())
plt.ylim(lat.min(),lat.max())
plt.colorbar()
plt.clim(0,8)
plt.title('Var_rand')

plt.show()
