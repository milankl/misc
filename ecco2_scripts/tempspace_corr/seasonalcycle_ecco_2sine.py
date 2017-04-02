## investigate the seasonal cycle in ECCO2

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc

#for fitting sine wave
from scipy.optimize import leastsq

## LOAD DATA

theta1 = np.load('python/ecco2/theta_sfc_NA.npy')
time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')

## LOCAL FUNCTIONS

def peval(x, p): #sine +  double freq sine
    return p[0] * np.sin(2 * np.pi * 1./365.25 * x + p[1]) \
    + p[2] * np.sin(2 * np.pi * 2./365.25 * x + p[3]) + p[4]

def residuals(p,y,x):
    return y - peval(x,p)

## LOOP OVER LON LAT
#preallocate
thvar = np.zeros(theta1.shape[1:])
thamp1 = np.zeros_like(thvar)
thamp2 = np.zeros_like(thvar)
thpha2 = np.zeros_like(thvar)
thpha2 = np.zeros_like(thvar)
thmean = np.zeros_like(thvar)

for ilon in range(lonna.shape[0]):
    print('procesing '+str(round(ilon*1.0/(lonna.shape[0])*100))+'%')
    for ilat in range(latna.shape[0]):
        if ~eccomask_NA[ilat,ilon]:
            #detrend
            theta1[:,ilat,ilon] = theta1[:,ilat,ilon] - np.polyval(np.polyfit(time,theta1[:,ilat,ilon],1),time)
            tstd = np.std(theta1[:,ilat,ilon])
            p0 = [tstd, 0, tstd/2., 0, np.mean(theta1[:,ilat,ilon])]
            thvar[ilat,ilon] = tstd
            
            plsq = leastsq(residuals, p0, args=(theta1[:,ilat,ilon], time))
            
            thamp1[ilat,ilon] = plsq[0][0]
            thpha1[ilat,ilon] = plsq[0][1]
            thamp2[ilat,ilon] = plsq[0][2]
            thpha2[ilat,ilon] = plsq[0][3]
            thmean[ilat,ilon] = plsq[0][4]
            
            if not(plsq[1] in [1,2,3,4]):
            	print('fitting not converging at '+str((ilat,ilon)))

            
np.save('python/ecco2/theta_seasonaldoublesine_NA_detr.npy',(thvar,thamp1,thpha1,thamp2,thpha2,thmean))
