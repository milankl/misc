## RECONSTRUCT T from dTdt?

import numpy as np
import matplotlib.pyplot as plt
import time as tictoc
from scipy.integrate import cumtrapz
exec(open('python/ecco2/local_functions.py').read())
exec(open('python/ecco2/colormap.py').read())

(eofs,pcs,eigs) = np.load('python/gyres/theta_eofs_lowres.npy')
(deofs,dpcs,deigs) = np.load('python/gyres/theta_eofs_highres.npy')

## TIME SCALES
(eofL, pctau) = np.load('python/gyres/eof_lowres_ltscales.npy')
(eofL2, pctau2) = np.load('python/gyres/eof_highres_ltscales.npy')

mmax = 200

fig, (ax1,ax2,ax3) = plt.subplots(3,sharex=True)

ax1.plot(pctau[:mmax])
ax1.plot(pctau2[:mmax])
ax1.set_ylabel('days')
ax1.set_title('Decorrelation time scale')

ax2.plot(eofL[:mmax])
ax2.plot(eofL2[:mmax])
ax2.set_ylabel('km')
ax2.set_title('Length scale')

ax3.plot(np.cumsum(eigs[:mmax])*100,label='EOF(Tlow)')
ax3.plot(np.cumsum(deigs[:mmax])*100,label='EOF(Thigh)')
ax3.set_ylabel('[%]')
ax3.set_xlabel('mode #')
ax3.legend(loc=4)

plt.show()

