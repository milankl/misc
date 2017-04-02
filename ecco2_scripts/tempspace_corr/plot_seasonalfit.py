## PLOT SEASONAL CYCLE ON MAP

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

## LOAD DATA

(thvar,thamp,thfreq,thpha,thmean) = np.load('python/ecco2/theta_seasonalsine_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

xx,yy = np.meshgrid(lonna,latna)

## APPLY MASK
thvar = np.ma.masked_array(thvar,mask=eccomask_NA)
thamp = np.ma.masked_array(thamp,mask=eccomask_NA)
thfreq = np.ma.masked_array(thfreq,mask=(eccomask_NA + (yy < 15.)))
thpha = np.ma.masked_array(thpha,mask=(eccomask_NA))
thmean = np.ma.masked_array(thmean,mask=eccomask_NA)

## display freq as period
thfreq[thfreq > 0] = 1/thfreq[thfreq > 0]
thamp = abs(thamp)
#thpha[thpha < 0] = thpha[thpha < 0] + 2*np.pi


## PLOTTING

plt.figure(1)
plt.pcolormesh(xx,yy,thvar)
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.colorbar()
plt.title('theta surface std')

plt.figure(2)
plt.pcolormesh(xx,yy,thamp/2/thvar)
plt.colorbar()
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.title('theta surface sine amp')

plt.figure(3)
plt.pcolormesh(xx,yy,thfreq)
plt.colorbar()
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.title('theta surface period')

plt.figure(4)
plt.pcolormesh(xx,yy,thpha)
plt.colorbar()
plt.clim((-.5,1))
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.title('theta surface phase')

plt.figure(5)
plt.pcolormesh(xx,yy,thmean)
plt.colorbar()
plt.xlim((lonna[0],lonna[-1]))
plt.ylim((latna[0],latna[-1]))
plt.title('theta surface mean')

"""
plt.figure(6)
m = Basemap(llcrnrlon=265,llcrnrlat=-5,urcrnrlon=345,urcrnrlat=55,projection='mill',\
rsphere=6371200.,resolution='c',area_thresh=5000)
m.drawcoastlines()
m.drawparallels(np.arange(0.,50.,5.),labels=[1,0,0,0],fontsize=10) 
m.drawmeridians(np.arange(265.,345.,5.),labels=[0,0,0,1],fontsize=10)
lons, lats = m(xx, yy)

m.pcolormesh(lons,lats,thmean)
plt.colorbar()
"""
plt.show()

