## VISUALIZE PATTERN

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['image.cmap'] = 'coolwarm'

(randfield2,mode) = np.load('python/ecco2/patterns/randpattern_4.npy')
(randfield1,mode) = np.load('python/ecco2/patterns/randpattern_1_ar4.npy')
theta = np.load('python/ecco2/theta_sfc_NA_ano_v2.npy')
eccomask_NA = np.load('python/ecco2/ecco_mask_NA.npy')

randfield2 = randfield2*2
randfield1 = randfield1*1.5

cmax = theta.std()*5


theta = np.ma.masked_array(theta,mask=np.array([eccomask_NA]*len(time)))

fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(16,5))

im1 = ax1.imshow(theta[0,:,:],interpolation='none',origin='lower')
im2 = ax2.imshow(randfield1[0,:,:],interpolation='none',origin='lower')
im3 = ax3.imshow(randfield2[0,:,:],interpolation='none',origin='lower')
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

ax2.set_title('SST random pattern, mode 1-100')
ax3.set_title('SST random pattern, mode 51-150')
ax1.set_title('SST ECCO2')

tmax = 600
timestamp = ax2.text(2,180,str(0)+' / '+str(tmax)+' days')

fig.canvas.draw()

plt.tight_layout()

for i in range(tmax):
    im1.set_array(theta[i,:,:])
    im2.set_array(randfield1[i,:,:])
    im3.set_array(randfield2[i,:,:])
    timestamp.set_text(str(i*3)+' / '+str(tmax*3)+' days')
    fig.canvas.draw()
    plt.tight_layout()
    plt.savefig('python/ecco2/frames/frame_%03d.png' % i, bbox_inches='tight')
    print('frame %03d saved to file.' % i)
    plt.pause(.001)
