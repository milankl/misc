## EOF ON THETA

import numpy as np
import matplotlib.pyplot as plt

## IMPORT DATA

time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
xx,yy = np.meshgrid(lonna,latna)

eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(thetaeofs,thetapcs) = np.load('python/ecco2/theta_eofs_1000.npy')

(eofL,pctau) = np.load('python/ecco2/eof_ltscales.npy')

# apply mask
thetaeofs = np.ma.masked_array(thetaeofs,mask=np.array([eccomask_NA]*1000))

## PLOTTING
mode = 9
years = 1992 + time/365.25

fig = plt.figure()

ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=2)
ax2 = plt.subplot2grid((3, 2), (2, 0), colspan=2)

im1 = ax1.pcolormesh(xx,yy,thetaeofs[mode-1,:,:],cmap='RdBu')
c = plt.colorbar(im1,ax=ax1)
cmax = abs(thetaeofs[mode-1,:,:]).max()
im1.set_clim(-cmax,cmax)
ax1.set_xlim(270,340)
ax1.set_ylim(0,50)
ax1.set_title('eof mode '+str(mode))

ax2.plot(years,thetapcs[:,mode])
ax2.set_title('Length scale '+str(eofL[mode-1])[:5]+'km, Time scale '+str(pctau[mode-1])[:5]+'days')

plt.show()
