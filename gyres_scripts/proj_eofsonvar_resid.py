## project EOFS onto Var(hi) - Var(lo)

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormap.py').read())
from matplotlib.colors import LogNorm

(eofs1,pcs1,eigs1) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs2,pcs2,eigs2) = np.load('python/gyres/theta_eofs_highres.npy')

tlo =  np.load('python/gyres/temp_lowres_sfc.npy')
thi =  np.load('python/gyres/temp_highres_sfc.npy')

(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

## DIFFERENCE IN VARIANCE

dvar = (thi.var(axis=0) - tlo.var(axis=0)).flatten()

n = 200
m = 200
s = 5

x = np.arange(0,n,s)
y = np.arange(1,m,s)
"""
q1 = np.zeros((n/s,m/s))
q2 = np.zeros((n/s,m/s))

for i,ii in zip(x,range(n/s)):
    print(i)
    for j,jj in zip(y,range(m/s)):
        e1 = (eofs1[i:i+j,...]**2).reshape((j,-1)).T
        e2 = (eofs2[i:i+j,...]**2).reshape((j,-1)).T
        
        q1[ii,jj] = np.linalg.norm(np.linalg.lstsq(e1,dvar)[1])
        q2[ii,jj] = np.linalg.norm(np.linalg.lstsq(e2,dvar)[1])

"""

v = np.arange(0,30000,1000)

fig,(ax1,ax2) = plt.subplots(1,2,sharey=True)

im1 = ax1.contourf(x,y,q1.T,v,extend='max')
plt.colorbar(im1,ax=ax1)
ax1.set_xlabel('mode_i')
ax1.set_ylabel('# modes')
ax1.set_title('norm(residuals) for T_low')


im2 = ax2.contourf(x,y,q2.T,v,extend='max')
plt.colorbar(im2,ax=ax2)
ax2.set_xlabel('mode_i')
ax2.set_title('norm(residuals) for T_high')


plt.show()

