## project EOFS onto Var(hi) - Var(lo)

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormap.py').read())

(eofs1,pcs1,eigs1) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs2,pcs2,eigs2) = np.load('python/gyres/theta_eofs_highres.npy')

tlo =  np.load('python/gyres/temp_lowres_sfc.npy')
thi =  np.load('python/gyres/temp_highres_sfc.npy')

(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

## DIFFERENCE IN VARIANCE

dvar = (thi.var(axis=0) - tlo.var(axis=0)).flatten()
mmax = 100

e1 = (eofs1[:mmax,...]**2).reshape((mmax,-1)).T
e2 = (eofs2[:mmax,...]**2).reshape((mmax,-1)).T

## 

q1 = np.linalg.lstsq(e1,dvar)[0]
q2 = np.linalg.lstsq(e2,dvar)[0]

fig,(ax1,ax2) = plt.subplots(2,sharex=True)
ax1.plot(abs(q1))
ax2.plot(abs(q2))
ax1.set_title('alpha_lo_i')
ax2.set_title('alpha_hi_i')

plt.figure(2)
plt.pcolormesh(thi.var(axis=0) - tlo.var(axis=0))
plt.colorbar()
plt.clim(dvar.min(),dvar.max())
plt.title('Var(T_hi) - Var(T_lo)')

plt.figure(3)
plt.pcolormesh(np.dot(e1,q1).reshape((128,128)))
plt.colorbar()
plt.clim(dvar.min(),dvar.max())
plt.title('Var(alpha*e)_lo')

plt.figure(4)
plt.pcolormesh(np.dot(e2,q2).reshape((128,128)))
plt.colorbar()
plt.clim(dvar.min(),dvar.max())
plt.title('Var(alpha*e)_high')

plt.show()
