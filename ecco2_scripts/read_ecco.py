## read ECCO2 only for one layer and one area

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncread
import time as tictoc

## local functions

def closest(vec,x):
    #returns the index of vector vec which contains a value that is closest to value x
    a = abs(vec - x)
    return np.where(a == np.min(a))[0][0]

## define paths and read file names

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

# area - 260 to 320 EAST, 30 to 50 NORTH
lonmin,lonmax = closest(lon,270),closest(lon,340)
latmin,latmax = closest(lat,0),closest(lat,50)

dep1 = 0 #only first level, i.e. 5m

## READ ALL NC FILES ONLY ONCE
# preallocate
theta1 = np.zeros((fnlist.shape[0],latmax-latmin+1,lonmax-lonmin+1))

tic0 = tictoc.time()
for ifn in range(len(fnlist)):
    dat = ncread(path+folder+fnlist[ifn])
    theta1[ifn,:,:] = dat.variables['THETA'][0,dep1,latmin:latmax+1,lonmin:lonmax+1]
    dat.close()
    print('reading '+str(np.round(float(ifn)/fnlist.shape[0]*100))+'%')

tic0 = tictoc.time() - tic0

np.save('python/ecco2/theta_sfc_NA.npy',theta1)
#theta1 = np.load('python/ecco2/theta_sfc_NA.npy')



