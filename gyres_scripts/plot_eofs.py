## EOF ON THETA

import numpy as np
import matplotlib.pyplot as plt

## IMPORT DATA

time,lat,lon = np.load('python/gyres/temp_lowres_dim.npy')
xx,yy = np.meshgrid(lon,lat)
time = time[::2]
time = (time - time[0])/24./3600.

(eofs1,pcs1,eigs1) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs2,pcs2,eigs2) = np.load('python/gyres/theta_eofs_highres.npy')
(eofL1,pctau1) = np.load('python/gyres/eof_lowres_ltscales.npy')
(eofL2,pctau2) = np.load('python/gyres/eof_highres_ltscales.npy')
years = time/365.25


## PLOTTING

mode = 1

fig1 = plt.figure(1)

ax11 = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=2)
ax12 = plt.subplot2grid((3, 2), (2, 0), colspan=2)

im11 = ax11.pcolormesh(xx,yy,eofs1[mode-1,:,:],cmap='RdBu')
c = plt.colorbar(im11,ax=ax11)
cmax = abs(eofs1[mode-1,:,:]).max()
im11.set_clim(-cmax,cmax)
ax11.set_xlim(lon.min(),lon.max())
ax11.set_ylim(lat.min(),lat.max())
ax11.set_title('eof mode '+str(mode))

ax12.plot(years,pcs1[:,mode])
ax12.set_title('Length scale '+str(eofL1[mode-1])[:5]+'km, Time scale '+str(pctau1[mode-1])[:5]+'days')
plt.tight_layout()

fig2 = plt.figure(2)

ax21 = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=2)
ax22 = plt.subplot2grid((3, 2), (2, 0), colspan=2)

im21 = ax21.pcolormesh(xx,yy,eofs2[mode-1,:,:],cmap='RdBu')
c = plt.colorbar(im21,ax=ax21)
cmax = abs(eofs2[mode-1,:,:]).max()
im21.set_clim(-cmax,cmax)
ax21.set_xlim(lon.min(),lon.max())
ax21.set_ylim(lat.min(),lat.max())
ax21.set_title('eof mode '+str(mode))

ax22.plot(years,pcs2[:,mode])
ax22.set_title('Length scale '+str(eofL2[mode-1])[:5]+'km, Time scale '+str(pctau2[mode-1])[:5]+'days')
plt.tight_layout()

plt.show()
