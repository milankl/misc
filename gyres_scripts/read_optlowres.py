## READ LOW RESOLUTION MITGCM RUN

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncread
import glob

path = '/home/jupiter/ocean/fenwick/baroclinic_gyre/baroclinic_128x128x4_uvTuvT_mk4/run_370/mnc_test_00'
#loop over 13,14,15,16

files = []

for folder in [13,14,15,16]:
    files.extend(glob.glob(path+str(folder)+'/state*.nc'))
    
x = np.zeros((len(files),33))
y = np.zeros((len(files),33))

#get dimensions

for i,fil in zip(range(len(files)),files):
    dat = ncread(fil)
    x[i,:] = dat.variables['X'][:]
    y[i,:] = dat.variables['Y'][:]
    
#read time only once
time = dat.variables['T'][:]
dat.close()
print('Dimensions read.')

x = np.sort(np.array(list(set([i for i in x.flatten()]))))
y = np.sort(np.array(list(set([i for i in y.flatten()]))))

## preallocate surface temperature
# use only 3000 last time steps (due to spin up)
tlength = 6000
temp = np.zeros((tlength,y.shape[0],x.shape[0]))

for i,fil in zip(range(len(files)),files):
    print('File no. '+str(i+1)+' of '+str(len(files)))
    dat = ncread(fil)
    datx = dat.variables['X'][:]
    daty = dat.variables['Y'][:]
    
    ys = np.where(y == daty[0])[0][0]
    xs = np.where(x == datx[0])[0][0]
    
    temp[:,ys:ys+daty.shape[0],xs:xs+datx.shape[0]]\
    = dat.variables['Temp'][:][-tlength:,0,:,:]

## get rid of boundaries
#adapt dimensions
# use only every other point in time to get daily resolution

bx = 128
by = 128

temp = temp[::2,:by,:bx]
time = time[-tlength::2]
time = (time - time[0])/24./3600.

y = y[:by]
x = x[:bx]

## SAVE

np.save('python/gyres/temp_optlowres_sfc.npy',temp)
#np.save('python/gyres/temp_lowres_dim.npy',(time,y,x))
print('Files written.')
