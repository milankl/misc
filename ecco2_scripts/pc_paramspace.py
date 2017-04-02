## ANALYSE AR PARAMETER SPACE FOR PRINCIPAL COMPONENTS

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/local_functions.py').read())

(thetaeofs,thetapcs) = np.load('python/ecco2/theta_eofs_1000.npy')

## AR2

n = 1000
alpha1 = np.zeros(n)
alpha2 = np.zeros(n)

for imode in range(n):
    (r1,r2) = acf(thetapcs[:,imode],3)[1:]

    alpha1[imode] = r1*(1-r2)/(1-r1**2)
    alpha2[imode] = (r2 - r1**2) / (1-r1**2)

## AR3
n = 1000
a1 = np.zeros(n)
a2 = np.zeros(n)
a3 = np.zeros(n)

for imode in range(n):
    (r1,r2,r3) = acf(thetapcs[:,imode],4)[1:]

    A = np.array([[1, r1, r2],[r1, 1, r1],[r2,r1,1]])
    (a1[imode],a2[imode],a3[imode]) = np.linalg.solve(A,np.array([r1,r2,r3]))


plt.figure(1)
plt.scatter(alpha1,alpha2,30,np.arange(n),alpha=.5)
plt.colorbar()
plt.plot([0,2,0],[-1,-1,1],'k')
plt.plot([0,1],[0,0],'k--')
plt.xlim(.4,2.1)
plt.ylim(-1.1,.1)
plt.xlabel('alpha1')
plt.ylabel('alpha2')

plt.figure(2)
plt.scatter(a1,a2,30+10*a3/np.std(a3),np.arange(n),alpha=.5)
plt.colorbar()
plt.plot([0,2,0],[-1,-1,1],'k')
plt.plot([0,1],[0,0],'k--')
plt.xlim(.4,2.1)
plt.ylim(-1.1,.1)
plt.xlabel('alpha1')
plt.ylabel('alpha2')
plt.show()
