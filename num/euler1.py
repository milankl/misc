## SOLVING ODEs WITH EULER FORWARD/MODIFIED

#ODE
#dy/dx = x/y, y(0) = 1

import numpy as np
import matplotlib.pyplot as plt

## analytic solution

xa = np.linspace(0,0.3,1000)
ya = np.sqrt(xa**2 + 1)

## solve with Euler forward

dx = 0.05
x = np.arange(0,0.3+dx,dx)


y1 = np.zeros_like(x)
y1[0] = 1 #initial condition

for i in range(len(x)-1):
    y1[i+1] = y1[i] + dx*x[i]/y1[i]

## solve with modified Euler

y2 = np.zeros_like(x)
y2[0] = 1 #initial condition

for i in range(len(x)-1):
    y2[i+1] = y2[i] + dx*x[i]/y2[i]
    y2[i+1] = y2[i] + .5*dx*(x[i]/y2[i] + x[i+1]/y2[i+1])

## plotting

plt.plot(xa,ya,'-',label='y analytic')
plt.plot(x,y,'s-',label='forw Euler')
plt.plot(x,y2,'s-',label='mod Euler')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc=2)
plt.title('Forward Euler vs. modified Euler')
plt.ylim(0.99,1.05)

plt.tight_layout()
plt.show()