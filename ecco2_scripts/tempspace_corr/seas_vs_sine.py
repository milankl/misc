## COMPARE SEASONAL CYCLE TO SINE WAVE

import numpy as np
import matplotlib.pyplot as plt
import datetime

#for fitting sine wave
from scipy.optimize import leastsq

def peval(x, p):
    return p[0] * np.sin(2 * np.pi * 1./365.25 * x + p[1]) \
    + p[2] * np.sin(2 * np.pi * 2./365.25 * x + p[3]) + p[4]

def residuals(p,y,x):
    return y - peval(x,p)

## LOAD DATA

theta = np.load('python/ecco2/theta_sfc_NA.npy')
time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')

## pick one point

lat1 = np.argmin(abs(latna - 48.))
lon1 = np.argmin(abs(lonna - 310.))
theta1 = theta[:,lat1,lon1]

## fit sine wave

p0 = [np.std(theta1), 0, np.std(theta1)/2., 0, np.mean(theta1)]
plsq = leastsq(residuals, p0, args=(theta1, time))
theta1s = peval(time,plsq[0])

## MULTI-YEAR 3DAILY AVERAGE

startdate = datetime.datetime(1992,1,1)
dt = [datetime.timedelta(d,0,0) for d in time]
dates = [startdate + d for d in dt]
years = [d.year for d in dates]
#include first and last index for looping later
newyearidx = np.hstack((0,np.where(np.diff(years))[0]+1))
newyeartype = [d.day for d in [dates[i] for i in newyearidx]]

#preallocate
ntype = [np.sum(np.array(newyeartype) == i) for i in [1,2,3]]
mycycle1 = np.zeros((ntype[0],np.diff(newyearidx)[3]))
mycycle2 = np.zeros((ntype[1],np.diff(newyearidx)[0]))
mycycle3 = np.zeros((ntype[2],np.diff(newyearidx)[2]))

#fill with nans for easier averaging as the last year is not complete
mycycle1.fill(np.nan)
mycycle2.fill(np.nan)
mycycle3.fill(np.nan)

i1,i2,i3 = 0,0,0
for ity in range(len(newyeartype)):
    
	if ity == (len(newyeartype) - 1): #to match correct indices for the last iteration
		lstart = newyearidx[ity]
		lend = None
	else:
		lstart,lend = newyearidx[ity],newyearidx[ity+1]
    
	if newyeartype[ity] == 1:
		mycycle1[i1,:] = theta1[lstart:lend]
		i1 = i1 + 1
	if newyeartype[ity] == 2:
		mycycle2[i2,:] = theta1[lstart:lend]
		i2 = i2 + 1
	if newyeartype[ity] == 3:
		mycycle3[i3,:theta1[lstart:lend].shape[0]] = theta1[lstart:lend]
		i3 = i3 + 1

# DO THE ACTUAL AVERAGE

mycycle1 = np.nanmean(mycycle1,axis=0)
mycycle2 = np.nanmean(mycycle2,axis=0)
mycycle3 = np.nanmean(mycycle3,axis=0)

# INTERWEAVE
# add 0  to have same size, kill afterwards

mycycle = np.vstack((mycycle1,mycycle2,np.hstack((mycycle3,0)))).flatten('F')[:-1]

## PLOT

time2 = np.arange(365)+1.5
sincycle = peval(time2,plsq[0])

plt.figure(1)
plt.plot(time2,mycycle,'r',time2,sincycle,'k')
plt.show()

