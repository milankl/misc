## PLOT TEMPORAL AUTOCORRELATIONS

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
exec(open('python/ecco2/colormap.py').read())

#tau1 = np.load('python/gyres/theta_lowres_acftau.npy')
tau1 = np.load('python/gyres/patterns/randpattern_ar5_20-120_high_acftau.npy')
tau2 = np.load('python/gyres/theta_highres_acftau.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

#correct lonna,latna for pcolor
dlon = lon[-1]-lon[-2]
dlat = lat[-1]-lat[-2]
lon = np.hstack((lon,lon[-1]+dlon))-dlon/2.
lat = np.hstack((lat,lat[-1]+dlat))-dlat/2.
xx,yy = np.meshgrid(lon,lat)

## PLOTTING

plt.figure(1,figsize=(10,6))
m = Basemap(llcrnrlon=-2,llcrnrlat=15,urcrnrlon=40,urcrnrlat=58,projection='mill',\
rsphere=6371200.,resolution='l',area_thresh=5000,fix_aspect=False)
m.drawparallels(np.arange(15.,65.,5.),labels=[1,0,0,0],fontsize=10) 
m.drawmeridians(np.arange(0.,45.,5.),labels=[0,0,0,1],fontsize=10)
lons, lats = m(xx, yy)

c = m.pcolormesh(lons,lats,tau1)
cc = plt.colorbar(c)
cc.set_label('days')
plt.clim((0,10))
plt.title('T_low: Decorrelation time scale')

plt.figure(2,figsize=(10,6))
m.drawparallels(np.arange(15.,65.,5.),labels=[1,0,0,0],fontsize=10) 
m.drawmeridians(np.arange(0.,45.,5.),labels=[0,0,0,1],fontsize=10)
c = m.pcolormesh(lons,lats,tau2)
cc = plt.colorbar(c)
cc.set_label('days')
plt.clim((0,10))
plt.title('T_high: Decorrelation time scale')

plt.show()
