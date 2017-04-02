## EOF ON THETA

import numpy as np
import matplotlib.pyplot as plt
from eofs.standard import Eof as eof
import time as tictoc
exec(open('python/ecco2/colormap.py').read())

## IMPORT DATA

theta =  np.load('python/gyres/temp_lowres_sfc.npy')
time,lat,lon = np.load('python/gyres/temp_lowres_dim.npy')

## EOF

#weights
wgts = np.sqrt(np.cos(np.deg2rad(lat)))[..., np.newaxis]
tic = tictoc.time()
thetasolv = eof(theta,weights=wgts)
print('EOFS of lowres data in '+str(tictoc.time() - tic)[:6]+'s.')

thetaeofs = thetasolv.eofs(neofs=1000,eofscaling=2)
thetapcs = thetasolv.pcs(npcs=1000,pcscaling=1)
thetaeigs = thetasolv.varianceFraction(neigs=1000)

np.save('python/gyres/theta_eofs_lowres.npy',(thetaeofs,thetapcs,thetaeigs))
print('Lowres EOFs written.')

## IMPORT DATA

theta =  np.load('python/gyres/temp_highres_sfc.npy')
time,lat,lon = np.load('python/gyres/temp_lowres_dim.npy')

## EOF

#weights
wgts = np.sqrt(np.cos(np.deg2rad(lat)))[..., np.newaxis]
tic = tictoc.time()
thetasolv = eof(theta,weights=wgts)
print('EOFS on highres data in '+str(tictoc.time() - tic)[:6]+'s.')

thetaeofs = thetasolv.eofs(neofs=1000,eofscaling=2)
thetapcs = thetasolv.pcs(npcs=1000,pcscaling=1)
thetaeigs = thetasolv.varianceFraction(neigs=1000)

np.save('python/gyres/theta_eofs_highres.npy',(thetaeofs,thetapcs,thetaeigs))
print('Highres EOFs written.')
