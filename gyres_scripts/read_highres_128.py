## READ LOW RESOLUTION MITGCM RUN

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncread
import glob

path = '/home/jupiter/ocean/fenwick/baroclinic_gyre/baroclinic_512x512x4/ensemble_arc/run_f3_trunc/'
ncfile = 'state.512_to_128.001.nc'

tstart = -3000
tend = None

dat = ncread(path+ncfile)
print('NC file found.')
x = dat.variables['X'][:128]
y = dat.variables['Y'][:128]
time = dat.variables['T'][:][tstart:tend]
temp = dat.variables['Temp'][:][tstart:tend,0,:128,:128]

print('Data read.')

## seconds to days
time = (time - time[0])/24./3600.

np.save('python/gyres/temp_highres_sfc.npy',temp)
np.save('python/gyres/temp_highres_dim.npy',(time,y,x))
print('Files written.')
