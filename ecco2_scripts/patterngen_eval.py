## EVALUATE PATTERNS BY SPATIAL AND TEMPORAL AUTOCORRELATION

## TEMPORAL AUTO CORRELATION

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
exec(open('python/ecco2/local_functions.py').read())

## LOAD DATA

randfield,mode = np.load('python/ecco2/patterns/randpattern_3.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

print('pattern_3 ...')

## AUTOCORRELATION LENGTH SCALE

actau = np.zeros(randfield.shape[1:])
lagmax = 100

for ilon in range(lonna.shape[0]):
    print('processing '+str(round(ilon*1.0/(lonna.shape[0])*100))+'%')
    for ilat in range(latna.shape[0]):
        if ~eccomask_NA[ilat,ilon]:
            #autocorrelation time scale (decrease to 1/e) in days
            try:
            	ac = findzc(acf(randfield[:,ilat,ilon],lagmax),1/np.exp(1))*3
            	if ac.shape[0] > 1:
            		print(str(ac)+'d crossings at '+str((lonna[ilon],latna[ilat])))
            	actau[ilat,ilon] = ac[0]
            except:
            	actau[ilat,ilon] = None
            	print('No crossings at '+str((lonna[ilon],latna[ilat])))

np.save('python/ecco2/patterns/randpattern_3_acftau.npy',actau)
