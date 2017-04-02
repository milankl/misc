## Energy Balance Model with forward Euler, modified Euler

import numpy as np
import matplotlib.pyplot as plt

## constants

global h,rho,C,a,S,e,s

h = 8.3e3 #height of troposphere
rho = 1.2 #density
C = 1e3 #specific heat capacity
a = 0.3 #albedo
S = 1367. #solar constant
e = 1 #emissivity
s = 5.67e-8 #stefan-boltzmann constant

## solve EBM numerically

def rhs(T):
    return ((1-a)*S/4. - e*s*T**4)/h/rho/C
    
def ebm(t):
    dt = t[1] - t[0]
    # preallocate 1: forw Euler, 2: mod. Euler
    T1,T2 = np.zeros_like(t),np.zeros_like(t)
    T1[0],T2[0] = 288,288 #initial condition
    
    for i in range(len(t)-1):
        #forward euler
        T1[i+1] = T1[i] + dt*rhs(T1[i])
        
        #modified euler
        T2[i+1] = T2[i] + dt*rhs(T2[i])
        T2[i+1] = T2[i] + .5*dt*(rhs(T2[i]) + rhs(T2[i+1]))
    
    #convert to celsius
    return T1-273.15,T2-273.15
    
#time stepping
td = np.linspace(0,200,30) #in days
t = td*3600*24 #in seconds (for solving the ODE)

T1,T2 = ebm(t)

e = .6 #change emissivity
T1e,T2e = ebm(t)

a = 0.26 #change albedo
T1a,T2a = ebm(t)

## plotting

plt.plot(td,T1,'gs-',label=r'forw Euler, $\epsilon = 1$, $\alpha = .3$')
plt.plot(td,T2,'g^-',label=r'mod Euler, $\epsilon = 1$, $\alpha = .3$')
plt.plot(td,T1e,'bs-',label=r'forw Euler, $\epsilon = .6$, $\alpha = .3$')
plt.plot(td,T2e,'b^-',label=r'mod Euler, $\epsilon = .6$, $\alpha = .3$')
plt.plot(td,T1a,'ks-',label=r'forw Euler, $\epsilon = .6$, $\alpha = .26$')
plt.plot(td,T2a,'k^-',label=r'mod Euler, $\epsilon = .6$, $\alpha = .26$')
plt.xlabel('time [days]')
plt.ylabel(r'temperature [$\degree{}$C]')
plt.legend(loc='center right')

plt.tight_layout()
plt.show()
    