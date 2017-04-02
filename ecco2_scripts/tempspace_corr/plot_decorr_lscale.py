## PLOT SPATIAL DECORRELATION

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
exec(open('python/ecco2/colormap.py').read())

s = 3

lscale,m = np.load('python/ecco2/patterns/decorr_lscale_'+str(s)+'.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

lscale = np.ma.masked_array(lscale,mask=m)
latna,lonna = latna[1:-1],lonna[1:-1]

#correct lonna,latna for pcolor
dlon = lonna[-1]-lonna[-2]
dlat = latna[-1]-latna[-2]
lonna = np.hstack((lonna,lonna[-1]+dlon))-dlon/2.
latna = np.hstack((latna,latna[-1]+dlat))-dlat/2.
xx,yy = np.meshgrid(lonna,latna)

## PLOTTING

plt.figure(1,figsize=(10,6))
m = Basemap(llcrnrlon=268,llcrnrlat=-2,urcrnrlon=342,urcrnrlat=52,projection='mill',\
rsphere=6371200.,resolution='l',area_thresh=5000,fix_aspect=False)
m.drawcoastlines()
m.drawparallels(np.arange(0.,50.,5.),labels=[1,0,0,0],fontsize=10) 
m.drawmeridians(np.arange(265.,345.,5.),labels=[0,0,0,1],fontsize=10)
lons, lats = m(xx, yy)

c = m.pcolormesh(lons,lats,lscale)
cc = plt.colorbar(c)
cc.set_label('km')
plt.clim(0,2000)
plt.title('Random pattern of Pot. Temp (EOF #31-60): Decorrelation length scale')
plt.show()


