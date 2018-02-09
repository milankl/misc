## logistic mapping


import numpy as np
import matplotlib.pyplot as plt

##

n = 20
N = 5
x1 = np.empty((N,n))
x1[:,0] = np.linspace(0.1,0.9,N)

x2 = np.empty((N,n))
x2[:,0] = np.linspace(0.1,0.9,N)

x3 = np.empty((N,n))
x3[:,0] = np.linspace(0.3,0.31,N)

r1 = 3
r2 = 3.449
r3 = 3.6

for i in range(n-1):
    x1[:,i+1] = r1 * x1[:,i] *( 1 - x1[:,i] )
    x2[:,i+1] = r2 * x2[:,i] *( 1 - x2[:,i] )
    x3[:,i+1] = r3 * x3[:,i] *( 1 - x3[:,i] )

plt.figure(1)
plt.subplot(3,1,1)

plt.plot(x1.T,'o--')
plt.grid()
plt.ylabel('x_n')

plt.subplot(3,1,2)

plt.plot(x2.T,'o--')
plt.grid()
plt.ylabel('x_n')

plt.subplot(3,1,3)

plt.plot(x3.T,'o--')
plt.grid()
plt.xlabel('n')
plt.ylabel('x_n')
plt.show()