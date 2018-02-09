# the Lorenz model is: (cyclical)
# dX[j]/dt=(X[j+1]-X[j-2])*X[j-1]-X[j]+F

import numpy as np
import matplotlib.pyplot as plt

## functions

def rk4(X):
    k1 = rhs(X)
    k2 = rhs(X+.5*h*k1)
    k3 = rhs(X+.5*h*k2)
    k4 = rhs(X+h*k3)
    return h/6. * (k1 + 2*k2 + 2*k3 + k4)
    
def rhs(X):
    return (X[(i+1)%J] - X[i-2]) * X[i-1] - X[i] + F


global J,i,h,F
J = 8 #the number of variables
i = np.arange(J)
h = 0.05 #the time step
F = 2. #the forcing

n = 1e3
t = np.linspace(0,(n-1)*h,n)
X = np.empty((J,n+1))

#the initial conditions (steady state)
X[:,0] = F
#perturbation
X[5,0] = .1

for j in range(int(n)):
  X[:,j+1] = X[:,j] + rk4(X[:,j])

## plotting

plt.plot(X[0,:])
plt.show()
