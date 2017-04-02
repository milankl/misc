## PLOT SEASONAL CYCLE ON MAP

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
exec(open('python/ecco2/colormaps.py').read())

## LOAD DATA

(thstd,thamp1,thpha1,thamp2,thpha2,thmean) = np.load('python/ecco2/theta_seasonaldoublesine_NA_detr.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

xx,yy = np.meshgrid(lonna,latna)

## APPLY MASK
thstd = np.ma.masked_array(thstd,mask=eccomask_NA)
thamp1 = np.ma.masked_array(thamp1,mask=eccomask_NA)
thamp2 = np.ma.masked_array(thamp2,mask=eccomask_NA)
thpha1 = np.ma.masked_array(thpha1,mask=eccomask_NA)
thpha2 = np.ma.masked_array(thpha2,mask=eccomask_NA)
thmean = np.ma.masked_array(thmean,mask=eccomask_NA)

## display freq as period
thamp1 = abs(thamp1)
thamp2 = abs(thamp2)
#thpha[thpha < 0] = thpha[thpha < 0] + 2*np.pi

thamp1r = (thamp1/np.sqrt(2))/thstd
thamp2r = (thamp2/np.sqrt(2))/thstd
thampr = np.sqrt((thamp1**2+thamp2**2)/2)/thstd

## PLOTTING

f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(16,10))

im1 = ax1.pcolormesh(xx,yy,thstd)
ax1.set_xlim((lonna[0],lonna[-1]))
ax1.set_ylim((latna[0],latna[-1]))
ax1.set_title('Pot. Temp: std')
plt.colorbar(im1,ax=ax1)
im1.set_clim(0,8)

im2 = ax2.pcolormesh(xx,yy,thamp1r)
ax2.set_xlim((lonna[0],lonna[-1]))
ax2.set_ylim((latna[0],latna[-1]))
ax2.set_title('Pot. Temp: sin1 std / std')
plt.colorbar(im2,ax=ax2)
im2.set_clim(0,1)

im3 = ax3.pcolormesh(xx,yy,thamp2r)
ax3.set_xlim((lonna[0],lonna[-1]))
ax3.set_ylim((latna[0],latna[-1]))
ax3.set_title('Pot. Temp: sine2 std / std')
plt.colorbar(im3,ax=ax3)
im3.set_clim(0,1)

im4 = ax4.pcolormesh(xx,yy,thampr)
ax4.set_xlim((lonna[0],lonna[-1]))
ax4.set_ylim((latna[0],latna[-1]))
ax4.set_title('Pot. Temp: sine1+2 std / std')
plt.colorbar(im4,ax=ax4)
im4.set_clim(0,1)

plt.show()










