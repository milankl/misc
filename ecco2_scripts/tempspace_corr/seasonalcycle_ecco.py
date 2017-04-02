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

def peval(x, p):
    return p[0] * np.sin(2 * np.pi * p[1] * x + p[2]) + p[3]

def residuals(p,y,x):
    return y - peval(x,p)

## LOOP OVER LON LAT
#preallocate
thvar = np.zeros(theta1.shape[1:])
thamp = np.zeros_like(thvar)
thpha = np.zeros_like(thvar)
thfreq = np.zeros_like(thvar)
thmean = np.zeros_like(thvar)

#initial guess for sine parameters
p0flag = True

for ilon in range(lonna.shape[0]):
    print('procesing '+str(round(ilon*1.0/(lonna.shape[0])*100))+'%')
    for ilat in range(latna.shape[0]):
        if ~eccomask_NA[ilat,ilon]:
            
            p0 = [np.std(theta1[:,ilat,ilon]), 1/365.25, 0, np.mean(theta1[:,ilat,ilon])]
            
            #processing local time series
            thvar[ilat,ilon] = np.std(theta1[:,ilat,ilon])
            
            tic1 = tictoc.time()
            plsq = leastsq(residuals, p0, args=(theta1[:,ilat,ilon], time))
            tic1 = tictoc.time() - tic1
            
            thamp[ilat,ilon] = plsq[0][0]
            thfreq[ilat,ilon] = plsq[0][1]
            thpha[ilat,ilon] = plsq[0][2]
            thmean[ilat,ilon] = plsq[0][3]
            
            if not(plsq[1] in [1,2,3,4]):
            	print('fitting not converging at '+str((ilat,ilon)))

            
np.save('python/ecco2/theta_seasonalsine_NA.npy',(thvar,thamp,thfreq,thpha,thmean))
