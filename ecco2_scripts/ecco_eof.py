## EOF ON THETA

import numpy as np
import matplotlib.pyplot as plt
from eofs.standard import Eof as eof
import time as tictoc
exec(open('python/ecco2/colormap.py').read())

## IMPORT DATA

theta =  np.load('python/ecco2/theta_sfc_NA_ano_v2.npy')
time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')

s = 1 # take only every s point in the beginning

xx,yy = np.meshgrid(lonna[::s],latna[::s])

#theta = np.ma.masked_array(theta,mask=np.array([eccomask_NA]*len(time)))
thetas = theta[:,::s,::s]

## EOF

#weights
wgts = np.sqrt(np.cos(np.deg2rad(latna[::s])))[..., np.newaxis]
tic = tictoc.time()
thetasolv = eof(thetas,weights=wgts)
print(tictoc.time() - tic)

thetaeofs = thetasolv.eofs(neofs=1000,eofscaling=2)
thetapcs = thetasolv.pcs(npcs=1000,pcscaling=1)
thetaeigs = thetasolv.varianceFraction(neigs=1000)

np.save('python/ecco2/theta_eofs_1000.npy',(thetaeofs,thetapcs,thetaeigs))
