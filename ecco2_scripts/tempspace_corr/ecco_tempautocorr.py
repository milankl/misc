## TEMPORAL AUTO CORRELATION

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

## LOCAL FUNCTIONS

def acf(x,l):
	return np.array([1]+[np.corrcoef(x[:-i],x[i:])[0,1] for i in range(1,l)])

def findzc(x,a):
	return np.where(abs(np.diff(np.ceil(x-a))) == 1)[0]

## LOAD DATA

theta1ano = np.load('python/ecco2/theta_sfc_NA_ano_v2.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

xx,yy = np.meshgrid(lonna,latna)

## AUTOCORRELATION LENGTH SCALE

thactau = np.zeros(theta1ano.shape[1:])
lagmax = 100

for ilon in range(lonna.shape[0]):
    print('processing '+str(round(ilon*1.0/(lonna.shape[0])*100))+'%')
    for ilat in range(latna.shape[0]):
        if ~eccomask_NA[ilat,ilon]:
            #autocorrelation time scale (decrease to 1/e) in days
            try:
            	ac = findzc(acf(theta1ano[:,ilat,ilon],lagmax),1/np.exp(1))*3
            	if ac.shape[0] > 1:
            		print(str(ac)+'d crossings at '+str((lonna[ilon],latna[ilat])))
            	thactau[ilat,ilon] = ac[0]
            except:
            	thactau[ilat,ilon] = None
            	print('No crossings at '+str((lonna[ilon],latna[ilat])))

np.save('python/ecco2/theta_acftau_detr.npy',thactau)
