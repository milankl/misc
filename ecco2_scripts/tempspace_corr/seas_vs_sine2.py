## COMPARE SEASONAL CYCLE TO SINE WAVE
# SECOND APPROACH

import numpy as np
import matplotlib.pyplot as plt
import datetime
import time as tictoc

#for fitting sine wave
from scipy.optimize import leastsq

def peval(x, p):
    return p[0] * np.sin(2 * np.pi * 1./365.25 * x + p[1]) \
    + p[2] * np.sin(2 * np.pi * 2./365.25 * x + p[3]) + p[4]

def residuals(p,y,x):
    return y - peval(x,p)

def rmean(x, b):
	a = (b-1)/2
	y = np.convolve(x,np.ones(b)/b)[a:-a]
	y[:a] = None
	y[-a:] = None
	return y

def rmeanp(x,b):
	return np.convolve(np.hstack((x[-(b-1)/2:],x,x[:(b-1)/2])),np.ones(b)/b)\
	[(b-1):-(b-1)]

def acf(x,l):
	return np.array([1]+[np.corrcoef(x[:-i],x[i:])[0,1] for i in range(1,l)])

## LOAD DATA

theta = np.load('python/ecco2/theta_sfc_NA.npy')
time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')

## pick one point

lat1 = np.argmin(abs(latna - 6.))
lon1 = np.argmin(abs(lonna - 335.))
theta1 = theta[:,lat1,lon1]

## detrend point

#theta1 = theta1 - np.polyval(np.polyfit(time,theta1,1),time)

## fit sine wave

p0 = [np.std(theta1), 0, np.std(theta1)/2., 0, np.mean(theta1)]
plsq = leastsq(residuals, p0, args=(theta1, time))
theta1s = peval(time,plsq[0])

timemod = np.round(time%365.25,2)
timey = np.linspace(0,365.,365)

## AVERAGE OVER SAME ENTRIES OF TIMEMOD

timeset = np.sort(np.array(list(set(timemod))))
cycle = np.array([np.mean(theta1[t == timemod]) for t in timeset])
cyclerm = rmeanp(cycle,11)

## REMOVE SEASONAL CYCLE - RUNNING MEAN

theta1ano = np.array([theta1[i] - cycle[timeset == timemod[i]] for i in range(theta1.shape[0])])[:,0]
theta1anorm = np.array([theta1[i] - cyclerm[timeset == timemod[i]] for i in range(theta1.shape[0])])[:,0]

print((np.corrcoef(theta1ano,theta1-theta1s)[1,0],\
np.corrcoef(theta1ano,theta1anorm)[1,0],\
np.corrcoef(theta1anorm,theta1-theta1s)[1,0]))

"""
plt.figure(1)
plt.plot(timemod,theta1,'.',markersize=2)
plt.plot(timey,peval(timey,plsq[0]),'r')
plt.plot(timeset[::3],cycle[::3],'g')
plt.plot(timeset,rmeanp(cycle,11),'k',linewidth=4)

plt.figure(2)
l = 2700
plt.plot(time[:l],acf(theta1-theta1s,l),label='sin rem')
plt.plot(time[:l],acf(theta1ano,l),label='cyc rem')
plt.plot(time[:l],acf(theta1anorm,l),label='rmean(cyc) rem')
plt.plot(time[:l],1/np.exp(1)*np.ones(l),'k--')
plt.plot(time[:l],-1/np.exp(1)*np.ones(l),'k--')
plt.grid()
plt.legend()

plt.figure(3)
plt.plot(time,theta1-theta1s,'r',label='sin rem')
plt.plot(time,theta1ano,'b',label='cyc rem')
plt.plot(time,theta1anorm,'g',label='rmean(cyc) rem')
plt.legend()
plt.show()
"""

f, (ax1,ax2) = plt.subplots(2,figsize=(10,10))
ax1.plot(timemod,theta1,'.',markersize=2)
ax1.plot(timey,peval(timey,plsq[0]),'r',linewidth=4)
ax1.plot(timeset,rmeanp(cycle,5),'k',linewidth=4)
ax1.set_xlim(0,365.25)
ax1.set_title('Seasonal cycle at '+str(latna[lat1])+'N, '+str(lonna[lon1])+'E')
ax1.set_ylabel(r'$\degree{} $C')
ax1.set_xlabel('day in year')

l = 300
ax2.plot(time[:l],acf(theta1-theta1s,l),'r',label='sine')
ax2.plot(time[:l],acf(theta1anorm,l),'k',label='rmean(cyc)')
ax2.plot(time[:l],acf(theta1,l),'b',label='unfiltered')
ax2.plot(time[:l],1/np.exp(1)*np.ones(l),'k--')
ax2.plot(time[:l],-1/np.exp(1)*np.ones(l),'k--')
ax2.set_xlabel('lag in days')
ax2.set_ylim(-0.6,1.)
ax2.legend()
ax2.set_title('Autocorrelation function')

plt.show()

