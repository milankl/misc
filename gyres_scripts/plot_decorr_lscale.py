## PLOT SPATIAL DECORRELATION

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
exec(open('python/ecco2/colormap.py').read())

#lscale1 = np.load('python/gyres/decorr_lscale_low_1.npy')
lscale1 = np.load('python/gyres/patterns/decorr_lscale_high_ar5_20-120_high.npy')
lscale2 = np.load('python/gyres/decorr_lscale_high_1.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

lat,lon = lat[1:-1],lon[1:-1]

#correct lon,lat for pcolor
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

c = m.pcolormesh(lons,lats,lscale1)
cc = plt.colorbar(c)
cc.set_label('km')
plt.clim((0,400))
plt.title('T_low: Decorrelation length scale')

plt.figure(2,figsize=(10,6))
m.drawparallels(np.arange(15.,65.,5.),labels=[1,0,0,0],fontsize=10) 
m.drawmeridians(np.arange(0.,45.,5.),labels=[0,0,0,1],fontsize=10)
c = m.pcolormesh(lons,lats,lscale2)
cc = plt.colorbar(c)
cc.set_label('km')
plt.clim((0,400))
plt.title('T_high: Decorrelation length scale')

plt.show()


