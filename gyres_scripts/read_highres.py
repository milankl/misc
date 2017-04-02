## READ LOW RESOLUTION MITGCM RUN

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncread
import glob

path = '/home/jupiter/ocean/fenwick/baroclinic_gyre/baroclinic_512x512x4/run5/mnc_test_0001/'
files = glob.glob(path+'state*.nc')

## as every dataset contains 94 time steps read only the last 32 files (skip the last as it is smaller)
# to obtain rouhgly 3000 time steps

nfiles = 32
tsteps = 94
nx = 512
ny = 512

#read only time first
files = files[-nfiles-1:-1]

#preallocate
time = np.zeros(nfiles*tsteps)
temp = np.zeros((nfiles*tsteps,ny,nx))

for i,fil in zip(range(len(files)),files):
    print('File no. '+str(i+1)+' of '+str(len(files)))
    dat = ncread(fil)
    time[i*tsteps:(i+1)*tsteps] = dat.variables['T'][:]
    temp[i*tsteps:(i+1)*tsteps,:,:] = dat.variables['Temp'][:][:,0,:ny,:nx]

#read grid only once
x = dat.variables['X'][:]
y = dat.variables['Y'][:]

## SAVE

time = (time - time[0])/24./3600.

np.save('python/gyres/temp_highres_sfc.npy',temp)
np.save('python/gyres/temp_highres_dim.npy',(time,y,x))
print('Files written.')
