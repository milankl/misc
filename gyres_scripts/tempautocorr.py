## TEMPORAL AUTO CORRELATION

import numpy as np
import time as tictoc

exec(open('python/ecco2/local_functions.py').read())

## LOAD DATA

#theta = np.load('python/gyres/temp_highres_sfc.npy')
theta = np.load('python/gyres/patterns/randpattern_ar5_20-120_high.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

## REMOVE MEAN
tic = tictoc.time()
theta = theta - theta.mean(axis=0)
print('Remove mean in '+str(tictoc.time()-tic)[:5]+'s.')

## AUTOCORRELATION LENGTH SCALE

thactau = np.zeros(theta.shape[1:])
lagmax = 200

tic = tictoc.time()
for ilon in range(lon.shape[0]):
    print('processing '+str(round(ilon*1.0/(lon.shape[0])*100))+'%')
    for ilat in range(lat.shape[0]):
        try:
        	ac = findzc(acfast(theta[:,ilat,ilon],lagmax),1/np.exp(1))
        	if ac.shape[0] > 1:
        		print(str(ac)+'d crossings at '+str((lon[ilon],lat[ilat])))
        	thactau[ilat,ilon] = ac[0]
        except:
        	thactau[ilat,ilon] = None
        	print('No crossings at '+str((lon[ilon],lat[ilat])))

print('Time scales in '+str(tictoc.time()-tic)[:5]+'s.')
np.save('python/gyres/patterns/randpattern_ar5_20-120_high_acftau.npy',thactau)
print('Files written.')
