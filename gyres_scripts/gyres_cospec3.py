## VARIANCE OF HIGH VS LOW

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormaps.py').read())
exec(open('python/ecco2/local_functions.py').read())
from matplotlib.colors import LogNorm
from scipy.signal import detrend

## load data

thi = np.load('python/gyres/temp_highres_sfc.npy')
tlo = np.load('python/gyres/temp_lowres_sfc.npy')
randfield = np.load('python/gyres/patterns/randpattern_ar5_0-100_high.npy')
randfield2 = np.load('python/gyres/patterns/randpattern_ar5_500-600_high.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

"""
## LOAD
(eofL1,pctau1) = np.load('python/gyres/eof_lowres_ltscales.npy')
(eofL2,pctau2) = np.load('python/gyres/eof_highres_ltscales.npy')

(eofs1,pcs1,eigs1) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs2,pcs2,eigs2) = np.load('python/gyres/theta_eofs_highres.npy')

tau1 = np.load('python/gyres/theta_lowres_acftau.npy')
tau2 = np.load('python/gyres/theta_highres_acftau.npy')
print('Data read.')
"""
"""
thi = detrend(thi,axis=0)
tlo = detrend(tlo,axis=0)
print('Data detrended.')
"""

dt = 1.
dy = (lat[1]-lat[0])*111.194
dx = (lon[1]-lon[0])*111.194*np.cos(2*np.pi*lat.mean()/360.)
"""
k,f,thihat = trispec(thi,dt,dy,dx)
print('FFT(T_high) done.')
k,f,tlohat = trispec(tlo,dt,dy,dx)
print('FFT(T_low) done.')
k,f,randhat = trispec(randfield,dt,dy,dx)
print('FFT(T_rand) done.')
k,f,rand2hat = trispec(randfield2,dt,dy,dx)
print('FFT(T_rand2) done.')

k = k[1:]
f = f[1:]
thihat = thihat[1:,1:]
tlohat = tlohat[1:,1:]
randhat = randhat[1:,1:]
rand2hat = rand2hat[1:,1:]
"""
v1 = 10.**np.arange(10)

fig1,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True,figsize=(14,10))
im1 = ax1.contourf(1/k,1/f,tlohat.T,v1,norm = LogNorm())
ax1.contour(1/k,1/f,tlohat.T,v1,colors='k',norm = LogNorm())
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim((1/k).min(),4e3)
ax1.set_ylim((1/f).min(),3e3)
ax1.set_title('Power spectra: SST low resolution')
ax1.set_ylabel('time scale [days]')

im2 = ax2.contourf(1/k,1/f,thihat.T,v1,norm = LogNorm())
ax2.contour(1/k,1/f,thihat.T,v1,colors='k',norm = LogNorm())
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim((1/k).min(),4e3)
ax2.set_ylim((1/f).min(),3e3)
ax2.set_title('SST high resolution')

im3 = ax3.contourf(1/k,1/f,randhat.T,v1,norm = LogNorm())
ax3.contour(1/k,1/f,randhat.T,v1,colors='k',norm = LogNorm())
plt.colorbar(im2,ax=(ax1,ax2,ax3,ax4))
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlim((1/k).min(),4e3)
ax3.set_ylim((1/f).min(),3e3)
ax3.set_title('SST random pattern, mode 1-100')
ax3.set_xlabel('length  scale [km]')
ax3.set_ylabel('time scale [days]')


im4 = ax4.contourf(1/k,1/f,rand2hat.T,v1,norm = LogNorm())
ax4.contour(1/k,1/f,rand2hat.T,v1,colors='k',norm = LogNorm())
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_xlim((1/k).min(),4e3)
ax4.set_ylim((1/f).min(),3e3)
ax4.set_title('SST random pattern, mode 500-600')
ax4.set_xlabel('length  scale [km]')

plt.show()



