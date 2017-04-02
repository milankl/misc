## EOF ON THETA

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc

## LOCAL FUNCTIONS
exec(open('python/ecco2/local_functions.py').read())

def nancorr(x,y):
    x = x - np.nanmean(x)
    y = y - np.nanmean(y)
    return np.nanmean(x*y)/np.std(x)/np.std(y)

## IMPORT DATA

time,lat,lon = np.load('python/gyres/temp_lowres_dim.npy')
xx,yy = np.meshgrid(lon,lat)

dl = lon[1]-lon[0]
dp = lat[1]-lat[0]
R = 6.371e6 #radius of earth

(eofs,pcs,eigs) = np.load('python/gyres/theta_eofs_highres.npy')

## ESTIMATE LENGTH SCALE
eofL = np.zeros(eofs.shape[0])

tic = tictoc.time()
for imode in range(1000):

    eord1 = eofs[imode,:,:].flatten('F')
    eord2 = eofs[imode,:,:].flatten('C')
    
    try:
        L1 = findzc(acfast(eord1,100),1/np.exp(1))[0]*dl*111.194
    except:
        L1 = np.nan
        print('No crossings for mode '+str(imode)+' in x.')
    
    try:
        L2 = findzc(acfast(eord2,100),1/np.exp(1))[0]*dp*111.194
    except:
        L2 = np.nan
        print('No crossings for mode '+str(imode)+' in y.')
    
    eofL[imode] = np.nanmean([L1,L2])

print('Length scales in '+str(tictoc.time() - tic)[:5]+'s.')
## TIME SCALE

pctau = np.zeros(eofs.shape[0])

tic = tictoc.time()
for imode in range(1000):
    try:
    	ac = findzc(acfast(pcs[:,imode],200),1/np.exp(1))
    	if ac.shape[0] > 1:
    		print(str(ac)+'d crossings for mode '+str(imode))
    	pctau[imode] = ac[0]
    except:
    	pctau[imode] = None
    	print('No crossings for mode '+str(imode))

print('Time scales in '+str(tictoc.time() - tic)[:5]+'s.')
np.save('python/gyres/eof_highres_ltscales.npy',(eofL,pctau))
print('LT scale file written.')


