## ECCO SPATIAL CORRELATION
## LOOP OVER LAT LON

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc

## LOAD DATA

#theta = np.load('python/gyres/temp_highres_sfc.npy')
theta = np.load('python/gyres/patterns/randpattern_ar5_20-120_high.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

## REMOVE MEAN
tic = tictoc.time()
theta = theta - theta.mean(axis=0)
theta = theta / theta.std(axis=0)
print('Normalise, standardise in '+str(tictoc.time()-tic)[:5]+'s.')

## SPATIAL CORRELATION FUNCTION
def mcorr(x,y):
    """correlation between vector x and 2d-array y. faster than np.corrcoef. x,y are assumed to have zero mean in time dimension (first dim) and variance = 1."""
    return (np.ma.dot(x,y) / (x.shape[0] - 1))

## CONSTANTS
mdeg = 111194.
acell = (lon[1]-lon[0])*(lat[1]-lat[0])*mdeg**2
w = np.array([acell*np.cos(np.pi / 180. * lat)]*lon.shape[0]).T

## PREALLOCATE

lscale = np.zeros(theta[0,1:-1,1:-1].shape)

for ilon in range(lon.shape[0])[1:-1]:
    tic = tictoc.time()
    print('processing '+str(round(ilon*1.0/(lon.shape[0])*100))+'%')
    for ilat in range(lat.shape[0])[1:-1]:
        #pic point
        theta1 = theta[:,ilat,ilon]
        
        ## INTELLIGENT SPATIAL EVALUATION
        # start around lat1,lon1
        
        #preallocate
        scorr = np.zeros(theta[0,:,:].shape)
        scorr[ilat,ilon] = 1.
        
        edges = np.array([True,True,True,True])   # up,right,down,left
        touch = np.array([False,False,False,False]) # boundary touched?
        rect = [ilon,ilon,ilat,ilat]    # lonmin,lonmax,latmin,latmax in index
        lim = 1/np.exp(1)
            
        while np.logical_and(edges,~touch).any():
            if edges[0] and not(touch[0]): #UP
                r = mcorr(theta1,theta[:,rect[3]+1,rect[0]-1:rect[1]+2])
                scorr[rect[3]+1,rect[0]-1:rect[1]+2] = r
                edges[0] = (r > lim).any()
                # reactivate left/right
                edges[3] = edges[3] or (r[0] > lim)
                edges[1] = edges[1] or (r[-1] > lim)
            
            if edges[1] and not(touch[1]): #RIGHT
                r = mcorr(theta1,theta[:,rect[2]-1:rect[3]+2,rect[1]+1])
                scorr[rect[2]-1:rect[3]+2,rect[1]+1] = r
                edges[1] = (r > lim).any()
                # reactivate up/down
                edges[2] = edges[2] or (r[0] > lim)
                edges[0] = edges[0] or (r[-1] > lim)        
            
            if edges[2] and not(touch[2]): #DOWN
                r = mcorr(theta1,theta[:,rect[2]-1,rect[0]-1:rect[1]+2])
                scorr[rect[2]-1,rect[0]-1:rect[1]+2] = r
                edges[2] = (r > lim).any()
                # reactivate left/right
                edges[3] = edges[3] or (r[0] > lim)
                edges[1] = edges[1] or (r[-1] > lim)
                   
            if edges[3] and not(touch[3]): #LEFT
                r = mcorr(theta1,theta[:,rect[2]-1:rect[3]+2,rect[0]-1])
                scorr[rect[2]-1:rect[3]+2,rect[0]-1] = r
                edges[3] = (r > lim).any()
                # reactivate up/down
                edges[2] = edges[2] or (r[0] > lim)
                edges[0] = edges[0] or (r[-1] > lim)
            
            #check boundaries
            touch[0] = ((rect[3] + 2) >= lat.shape[0])
            touch[1] = ((rect[1] + 2) >= lon.shape[0])
            touch[2] = (rect[2] < 2)
            touch[3] = (rect[0] < 2)
            
            if edges[0] and not(touch[0]): rect[3] = rect[3] + 1
            if edges[1] and not(touch[1]): rect[1] = rect[1] + 1
            if edges[2] and not(touch[2]): rect[2] = rect[2] - 1
            if edges[3] and not(touch[3]): rect[0] = rect[0] - 1

        ## CALCULATE DECORRELATION LENGTH SCALE [KM]
        lscale[(ilat-1),(ilon-1)] = np.sqrt(np.sum((scorr > lim)*w)/np.pi) / 1e3
        
    print "%.3f" % (tictoc.time() - tic)


np.save('python/gyres/patterns/decorr_lscale_high_ar5_20-120_high.npy',lscale)



