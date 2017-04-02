## EOF ON THETA

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc

## LOCAL FUNCTIONS

def acf(x,l):
	return np.array([1]+[nancorr(x[:-i],x[i:]) for i in range(1,l)])

def findzc(x,a):
	return np.where(abs(np.diff(np.ceil(x-a))) == 1)[0]

def nancorr(x,y):
    x = x - np.nanmean(x)
    y = y - np.nanmean(y)
    return np.nanmean(x*y)/np.std(x)/np.std(y)

## IMPORT DATA

time,latna,lonna = np.load('python/ecco2/ecco_dim_NA.npy')
xx,yy = np.meshgrid(lonna,latna)

dl = lonna[1]-lonna[0]
dp = latna[1]-latna[0]

R = 6.371e6

eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(eofs,pcs) = np.load('python/ecco2/theta_eofs_1000.npy')

# apply mask
eofs = np.ma.masked_array(eofs,mask=np.array([eccomask_NA]*1000))

## weights
wgts = np.cos(np.deg2rad(yy))
swgts = np.sum(wgts)

## AMPLITUDE

amp = (abs(eofs)).max(axis=2).max(axis=1)

## ESTIMATE LENGTH SCALE
eofL = np.zeros(eofs.shape[0])

tic = tictoc.time()
for imode in range(1000):

    eord1 = eofs[imode,:,:].flatten('F')
    eord2 = eofs[imode,:,:].flatten('C')
    
    try:
        L1 = findzc(acf(eord1,int(50*np.exp(-.03*imode) +10)),1/np.exp(1))[0]*2*dl*111.194
    except:
        L1 = np.nan
    
    try:
        L2 = findzc(acf(eord2,int(50*np.exp(-.03*imode) +10)),1/np.exp(1))[0]*2*dp*111.194
    except:
        L2 = np.nan
    
    eofL[imode] = np.nanmean([L1,L2])

print(tictoc.time() - tic)
## TIME SCALE

pctau = np.zeros(eofs.shape[0])

tic = tictoc.time()
for imode in range(1000):
    try:
    	ac = findzc(acf(pcs[:,imode],int(100*np.exp(-.02*imode) +10)),1/np.exp(1))*3
    	if ac.shape[0] > 1:
    		print(str(ac)+'d crossings for mode '+str(imode))
    	pctau[imode] = ac[0]
    except:
    	pctau[imode] = None
    	print('No crossings for mode '+str(imode))

print(tictoc.time() - tic)

np.save('python/ecco2/eof_ltscales.npy',(eofL,pctau))

"""
## PLOTTING STATISTICS

fig, (ax1,ax2, ax3) = plt.subplots(3,sharex=True)

ax1.plot(amp)
ax1.set_ylabel(r'$\degree{}$C')
ax1.set_title('Amplitude')

ax2.plot(pctau)
ax2.set_ylabel('days')
ax2.set_title('Decorrelation time scale')

ax3.plot(eofL)
ax3.set_ylabel('km')
ax3.set_title('Length scale')

plt.show()
"""

