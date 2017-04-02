## read ECCO2 dimensions

import numpy as np
from netCDF4 import Dataset as ncread

## local functions

def closest(vec,x):
    #returns the index of vector vec which contains a value that is closest to value x
    return np.argmin(abs(vec - x))

path = '/network/aopp/chaos/pred/cooper/data/ecco2.jpl.nasa.gov/data1/cube/cube92/lat_lon/quart_90S_90N/'
folder = 'THETA.nc/'
fnlist = open(path+folder+'fileNames.txt','r').readlines()
#remove the \n from every string
fnlist = [fn[:-1] for fn in fnlist]
#remove the fileNames.txt from that list
fnlist = np.delete(fnlist,np.where(np.array(fnlist) == 'fileNames.txt')[0][0])

## create a mask and get dimensions

dat1 = ncread(path+folder+fnlist[0])
eccomask = dat1.variables['THETA'][:].mask.squeeze()
lat = dat1.variables['LATITUDE_T'][:]
lon = dat1.variables['LONGITUDE_T'][:]
dep = dat1.variables['DEPTH_T'][:]
dat1.close()

#time in days since 1992-01-01
time = np.arange(len(fnlist))*3 + 1.5

## choose area - NORTH ATLANTIC
lonmin,lonmax = closest(lon,270),closest(lon,340)
latmin,latmax = closest(lat,0),closest(lat,50)

#for NA
latna = lat[latmin:latmax+1]
lonna = lon[lonmin:lonmax+1]

#mask for NA
eccomask_NA = eccomask[0,latmin:latmax+1,lonmin:lonmax+1]

np.save('python/ecco2/ecco_mask.npy',eccomask)
np.save('python/ecco2/ecco_mask_NA.npy',eccomask_NA)
np.save('python/ecco2/ecco_dim.npy',(time,dep,lat,lon))
np.save('python/ecco2/ecco_dim_NA.npy',(time,latna,lonna))
