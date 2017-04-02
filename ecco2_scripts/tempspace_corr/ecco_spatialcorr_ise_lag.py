## ECCO SPATIAL CORRELATION
## INCLUDING THE LAG

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc

## LOAD DATA

#theta = np.load('python/ecco2/theta_sfc_NA_ano_detr.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')
(time,latna,lonna) = np.load('python/ecco2/ecco_dim_NA.npy')

## SPATIAL CORRELATION FUNCTION
def mcorr(x,y):
    """correlation between vector x and 2d-array y. faster than np.corrcoef. x,y are assumed to have zero mean in time dimension (first dim)."""
    return ((np.ma.dot(x,y) / (x.shape[0] - 1) / y.std(axis=0)) / x.std())

## REMOVE MEAN
theta = theta - theta.mean(axis=0)

## define lags
lagmax = 10
lags = np.hstack((np.arange(lagmax+1),-np.arange(1,lagmax+1)))

## pic point
lat1 = np.argmin(abs(latna - 32.))
lon1 = np.argmin(abs(lonna - 283.))
theta1 = theta[:,lat1,lon1]

##preallocate
scorrmap = np.ma.masked_array(np.zeros(theta[:(2*lagmax+1),:,:].shape),\
mask=np.array([eccomask_NA]*(2*lagmax+1)))

# centers where to start ISE
#cent = [[lat1,lon1]]*(2*lagmax+1)
cent = [[]]*(2*lagmax+1)

for ilag in lags:
    print('processing '+str(round(ilag*1.0/(2*lagmax+1)*100))+'%')            
    
    ## INTELLIGENT SPATIAL EVALUATION
    tic = tictoc.time()
    
    ## APPLY LAG
    if ilag > 0:
        theta1 = theta[:-ilag,lat1,lon1]
        thetal = theta[ilag:,:,:]
    elif ilag < 0:
        theta1 = theta[-ilag:,lat1,lon1]
        thetal = theta[:ilag,:,:]
    else:
        theta1 = theta[:,lat1,lon1]
        thetal = theta
    
    ## NORMALIZE
    theta1 = theta1 - theta1.mean()
    thetal = thetal - thetal.mean(axis=0)
    
    lim = .15
    edges = np.array([True,True,True,True])   # up,right,down,left
    touch = np.array([False,False,False,False]) # boundary touched?
    
    # lonmin,lonmax,latmin,latmax in index (start values)
    if (ilag == 0) or (ilag == -1):
        rect = [lon1,lon1,lat1,lat1]
    elif ilag > 0:
        rect = [cent[ilag-1][1],cent[ilag-1][1],cent[ilag-1][0],cent[ilag-1][0]]
    elif ilag < -1:
        rect = [cent[ilag+1][1],cent[ilag+1][1],cent[ilag+1][0],cent[ilag+1][0]]
        
    #print(rect)
    # CALCULATE CENTERED CORRELATION
    scorrmap[ilag,rect[2],rect[0]] =\
    mcorr(theta1,thetal[:,rect[2],rect[0]])
        
    while np.logical_and(edges,~touch).any():
        if edges[0] and not(touch[0]): #UP
            r = mcorr(theta1,thetal[:,rect[3]+1,rect[0]-1:rect[1]+2])
            scorrmap[ilag,rect[3]+1,rect[0]-1:rect[1]+2] = r
            edges[0] = (r > lim).any()
            # reactivate left/right
            edges[3] = edges[3] or (r[0] > lim)
            edges[1] = edges[1] or (r[-1] > lim)
        
        if edges[1] and not(touch[1]): #RIGHT
            r = mcorr(theta1,thetal[:,rect[2]-1:rect[3]+2,rect[1]+1])
            scorrmap[ilag,rect[2]-1:rect[3]+2,rect[1]+1] = r
            edges[1] = (r > lim).any()
            # reactivate up/down
            edges[2] = edges[2] or (r[0] > lim)
            edges[0] = edges[0] or (r[-1] > lim)        
        
        if edges[2] and not(touch[2]): #DOWN
            r = mcorr(theta1,thetal[:,rect[2]-1,rect[0]-1:rect[1]+2])
            scorrmap[ilag,rect[2]-1,rect[0]-1:rect[1]+2] = r
            edges[2] = (r > lim).any()
            # reactivate left/right
            edges[3] = edges[3] or (r[0] > lim)
            edges[1] = edges[1] or (r[-1] > lim)
               
        if edges[3] and not(touch[3]): #LEFT
            r = mcorr(theta1,thetal[:,rect[2]-1:rect[3]+2,rect[0]-1])
            scorrmap[ilag,rect[2]-1:rect[3]+2,rect[0]-1] = r
            edges[3] = (r > lim).any()
            # reactivate up/down
            edges[2] = edges[2] or (r[0] > lim)
            edges[0] = edges[0] or (r[-1] > lim)
        
        #check boundaries
        touch[0] = ((rect[3] + 2) >= latna.shape[0])
        touch[1] = ((rect[1] + 2) >= lonna.shape[0])
        touch[2] = (rect[2] < 2)
        touch[3] = (rect[0] < 2)
        
        if edges[0] and not(touch[0]): rect[3] = rect[3] + 1
        if edges[1] and not(touch[1]): rect[1] = rect[1] + 1
        if edges[2] and not(touch[2]): rect[2] = rect[2] - 1
        if edges[3] and not(touch[3]): rect[0] = rect[0] - 1

    ## SET NEW STARTING POINT
    #if (ilag > 0) and (ilag < lagmax):
    #    cent[ilag+1] = np.unravel_index(np.argmax(scorrmap[ilag,:,:]),scorrmap.shape[1:])
    #if (ilag < 0) and (ilag > -lagmax):
    #    cent[ilag-1] = np.unravel_index(np.argmax(scorrmap[ilag,:,:]),scorrmap.shape[1:])
    
    cent[ilag] = np.unravel_index(np.argmax(scorrmap[ilag,:,:]),scorrmap.shape[1:])
    
    print "%.3f" % (tictoc.time() - tic)
    
## PLOTTING

cents = np.array(cent)
cents = np.vstack((cents[-lagmax:,:],cents[:lagmax+1,:]))
centtext = [str(i) for i in np.arange(-lagmax,lagmax+1)]

#correct lonna,latna for pcolor
dlon = lonna[-1]-lonna[-2]
dlat = latna[-1]-latna[-2]
lonna2 = np.hstack((lonna,lonna[-1]+dlon))-dlon/2.
latna2 = np.hstack((latna,latna[-1]+dlat))-dlat/2.
xx,yy = np.meshgrid(lonna2,latna2)

plt.figure(1)
plt.pcolormesh(xx,yy,scorrmap[-lagmax,:,:])
plt.colorbar()
plt.clim(0,1)

plt.figure(2)
plt.pcolormesh(xx,yy,scorrmap[0,:,:])
plt.plot(lonna[cents[:,1]],latna[cents[:,0]],'w--')
[plt.text(lonna[cents[i,1]],latna[cents[i,0]],centtext[i]) for i in range(2*lagmax+1)]
plt.colorbar()
plt.clim(0,1)

plt.figure(3)
plt.pcolormesh(xx,yy,scorrmap[lagmax,:,:])
plt.colorbar()
plt.clim(0,1)

plt.show()


