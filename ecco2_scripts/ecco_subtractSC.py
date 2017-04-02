## SUBTRACT SEASONAL CYCLE AS DOUBLE SINE

import numpy as np

#for fitting sine wave
from scipy.optimize import leastsq

## LOAD DATA

theta = np.load('python/ecco2/theta_sfc_NA.npy')
time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')

## LOCAL FUNCTIONS

## UPDATE: choose statistical model as sine, sine2, sine3, sine4, trend, constant

def peval(x, p):
    return p[0] * np.sin(2 * np.pi * 1./365.25 * x + p[1]) \
    + p[2] * np.sin(2 * np.pi * 2./365.25 * x + p[3]) \
    + p[4] * np.sin(2 * np.pi * 3./365.25 * x + p[5]) \
    + p[6] * np.sin(2 * np.pi * 4./365.25 * x + p[7]) \
    + p[8]*x + p[9]

def residuals(p,y,x):
    return y - peval(x,p)

## LOOP OVER LON LAT
#preallocate
nparam = 10 #number of parameters to fit
theta_ano = np.zeros_like(theta)
stat_param = np.zeros_like(theta[:10,:,:])

for ilon in range(lonna.shape[0]):
    print('processing '+str(round(ilon*1.0/(lonna.shape[0])*100))+'%')
    for ilat in range(latna.shape[0]):
        if ~eccomask_NA[ilat,ilon]:
            
            #detrend and subtract mean
            theta1 = theta[:,ilat,ilon]
            
            tstd = np.std(theta1)
            p0 = [tstd, 0, tstd/2., 0, tstd/3., 0, tstd/4., 0, 0, np.mean(theta1)]
                        
            plsq = leastsq(residuals, p0, args=(theta1, time))
            theta_ano[:,ilat,ilon] = theta1 - peval(time,plsq[0])
            stat_param[:,ilat,ilon] = plsq[0]
            
            if not(plsq[1] in [1,2,3,4]):
            	print('fitting not converging at '+str((ilat,ilon)))

            
np.save('python/ecco2/theta_sfc_NA_ano_v2.npy',theta_ano)
np.save('python/ecco2/theta_stat_param_v2.npy',stat_param)
