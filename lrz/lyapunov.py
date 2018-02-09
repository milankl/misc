#LYAPUNOV FRACTAL
import numpy as np
import matplotlib.pyplot as plt

n = 20
xl = 200
yl = 200

a = np.linspace(0.5,4,xl)
b = np.linspace(0.5,4,yl)

aa,bb = np.meshgrid(a,b)
S = 'AABAB'
l = len(S)

rn = np.empty(n)
x = np.empty(n)
lam = np.empty((xl,yl))

for i in range(xl):
    for ii in range(yl):
        xpos = aa[i,ii]
        ypos = bb[i,ii]
        
        x[0] = 0.5
        
        for iii in range(n-1):
            if S[np.mod(iii,l)] == 'A':
                rn[iii] = xpos
            else:
                rn[iii] = ypos
            
            x[iii+1] = rn[iii]*x[iii]*(1-x[iii])
        
        if S[np.mod(n-1,l)] == 'A':
            rn[n-1] = xpos
        else:
            rn[n-1] = ypos
        
        
        #now calculate lyapunov lambda
        lam[i,ii] = 1/float(n)*sum(np.log(abs(rn[1:]*(1-2*x[1:]))))

##plot

cmappi=plt.get_cmap('ocean')

fig, ax = plt.subplots()
cc = ax.contourf(aa,bb,lam,np.linspace(-2,.5,100),cmap=cmappi)
fig.colorbar(cc)
plt.show()