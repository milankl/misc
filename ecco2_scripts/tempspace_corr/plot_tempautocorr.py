## PLOT TEMPORAL AUTOCORRELATIONS

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
exec(open('python/ecco2/colormap.py').read())

tau = np.load('python/ecco2/theta_acftau_detr.npy')
#tau = np.load('python/ecco2/patterns/randpattern_1_acftau.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

#correct lonna,latna for pcolor
dlon = lonna[-1]-lonna[-2]
dlat = latna[-1]-latna[-2]
lonna = np.hstack((lonna,lonna[-1]+dlon))-dlon/2.
latna = np.hstack((latna,latna[-1]+dlat))-dlat/2.
xx,yy = np.meshgrid(lonna,latna)

#APPLY MASK
tau = np.ma.masked_array(tau,mask=(eccomask_NA + np.isnan(tau)))

## PLOTTING

plt.figure(1,figsize=(10,6))
m = Basemap(llcrnrlon=268,llcrnrlat=-2,urcrnrlon=342,urcrnrlat=52,projection='mill',\
rsphere=6371200.,resolution='l',area_thresh=5000,fix_aspect=False)
m.drawcoastlines()
m.drawparallels(np.arange(0.,50.,5.),labels=[1,0,0,0],fontsize=10) 
m.drawmeridians(np.arange(265.,345.,5.),labels=[0,0,0,1],fontsize=10)
lons, lats = m(xx, yy)

c = m.pcolormesh(lons,lats,tau)
cc = plt.colorbar(c)
cc.set_label('days')
plt.clim((0,180))
#plt.title('Random pattern of Pot. Temp (EOF #1-100): Decorrelation time scale')
plt.title('Pot. Temp: Decorrelation time scale')
plt.show()
