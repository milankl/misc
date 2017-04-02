## VISUALIZE PATTERN

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['image.cmap'] = 'RdBu'

randfield = np.load('python/gyres/patterns/randpattern_ar5_0-100_high.npy')
thi =  np.load('python/gyres/temp_highres_sfc.npy')
tlo =  np.load('python/gyres/temp_lowres_sfc.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

thi = thi - thi.mean(axis=0)
tlo = (tlo - tlo.mean(axis=0))*1.5

cmax = randfield.std()*4

fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,5))
im1 = ax1.imshow(randfield[0,:,:],interpolation='none',origin='lower')
im2 = ax2.imshow(tlo[0,:,:],interpolation='none',origin='lower')
im3 = ax3.imshow(thi[0,:,:],interpolation='none',origin='lower')
im1.set_clim(-cmax,cmax)
im2.set_clim(-cmax,cmax)
im3.set_clim(-cmax,cmax)
#plt.colorbar(im1,ax=(ax1,ax2,ax3),orientation='vertical',aspect=40)

ax1.set_ylabel('latitude')
ax1.set_xlabel('longitude')
ax2.set_xlabel('longitude')
ax3.set_xlabel('longitude')
ax1.set_xticks([])
ax2.set_xticks([])
ax3.set_xticks([])

ax1.set_yticks([])
ax2.set_yticks([])
ax3.set_yticks([])

ax1.set_title('SST random pattern, mode 1-100')
ax2.set_title('SST low resolution')
ax3.set_title('SST high resolution')

timestamp = ax2.text(5,120,str(0)+' / '+str(tmax)+' days')

fig.canvas.draw()

plt.tight_layout()

tmax = 600

for i in range(tmax):
    im1.set_array(randfield[i,:,:])
    im2.set_array(tlo[i,:,:])
    im3.set_array(thi[i,:,:])
    timestamp.set_text(str(i)+' / '+str(tmax)+' days')
    fig.canvas.draw()
    plt.savefig('python/gyres/frames/frame_%03d.png' % i, bbox_inches='tight')
    print('frame %03d saved to file.' % i)
    plt.pause(.001)
