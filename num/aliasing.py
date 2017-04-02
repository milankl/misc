## NUMERICS 3 - ALIASING ERROR

import numpy as np
import matplotlib.pyplot as plt

## grid

dx = 1
x = np.arange(0,10*dx,dx)
xa = np.linspace(0,max(x),1000)

kmax = np.pi/dx

k1 = 1
k2 = 2*kmax - k1

y1 = np.cos(k1*x)
y2 = np.cos(k2*x)

y1a = np.cos(k1*xa)
y2a = np.cos(k2*xa)


## plotting


plt.plot(xa,y1a,label=r'$k = 1, dx \to 0$')
plt.plot(xa,y2a,label=r'$k = 2\pi-1, dx \to 0$')
plt.plot(x,y1,'ko-',lw=2,label=r'$k = 1, dx = 1$')
plt.plot(x,y2,'rd--',lw=1,label=r'$k = 2\pi - 1, dx = 1$')
plt.title('Aliasing error')
plt.xlabel('x')
plt.ylabel(r'$\cos(kx)$')
plt.legend(loc='best')
plt.tight_layout()
plt.show()