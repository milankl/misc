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
topt = np.load('python/gyres/temp_optlowres_sfc.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

dt = 1.
dy = (lat[1]-lat[0])*111.194
dx = (lon[1]-lon[0])*111.194*np.cos(2*np.pi*lat.mean()/360.)

k,f,thihat = trispec(thi,dt,dy,dx)
print('FFT(T_high) done.')
k,f,tlohat = trispec(tlo,dt,dy,dx)
print('FFT(T_low) done.')
k,f,topthat = trispec(topt,dt,dy,dx)
print('FFT(T_optlow) done.')

k = k[1:]
f = f[1:]
thihat = thihat[1:,1:]
tlohat = tlohat[1:,1:]
topthat = topthat[1:,1:]

v1 = 10.**np.arange(10)

fig1,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
im1 = ax1.contourf(1/k,1/f,tlohat.T,v1,norm = LogNorm())
ax1.contour(1/k,1/f,tlohat.T,v1,colors='k',norm = LogNorm())
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim((1/k).min(),4e3)
ax1.set_ylim((1/f).min(),3e3)
ax1.set_title('T_low')
ax1.set_ylabel('time scale [days]')
ax1.set_xlabel('length scale [km]')

im2 = ax2.contourf(1/k,1/f,thihat.T,v1,norm = LogNorm())
ax2.contour(1/k,1/f,thihat.T,v1,colors='k',norm = LogNorm())

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim((1/k).min(),4e3)
ax2.set_ylim((1/f).min(),3e3)
ax2.set_title('T_high')
ax2.set_xlabel('length scale [km]')

im3 = ax3.contourf(1/k,1/f,topthat.T,v1,norm = LogNorm())
ax3.contour(1/k,1/f,topthat.T,v1,colors='k',norm = LogNorm())
plt.colorbar(im3,ax=(ax1,ax2,ax3))
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlim((1/k).min(),4e3)
ax3.set_ylim((1/f).min(),3e3)
ax3.set_title('T_opt')
ax3.set_xlabel('length  scale [km]')

v2 = 10.**np.linspace(-2,2,10)

im4 = ax4.contourf(1/k,1/f,topthat.T/thihat.T,v2,norm = LogNorm())
ax4.contour(1/k,1/f,topthat.T/thihat.T,v2,colors='k',norm = LogNorm())
plt.colorbar(im4,ax=(ax4))
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_xlim((1/k).min(),4e3)
ax4.set_ylim((1/f).min(),3e3)
ax4.set_title('T_opt/T_high')
ax4.set_xlabel('length  scale [km]')

plt.show()



