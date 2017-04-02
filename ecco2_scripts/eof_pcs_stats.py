## PRINCIPAL COMPENENTS VS WHITE NOISE

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormap.py'))
from scipy.stats import norm


(eofs,pcs,eigs) = np.load('python/ecco2/theta_eofs_1000.npy')

s = stats.skew(pcs,axis=0)
k = stats.kurtosis(pcs,axis=0)

#HISTOGRAM

bins = np.linspace(-6,6,150)
dbin = bins[1]-bins[0]
x = np.hstack((bins.min()-dbin,bins+dbin))

gauss = np.exp(-.5*x**2)/np.sqrt(2*np.pi)
mhist1 = np.histogram(pcs[:,8],bins)[0]
mhist2 = np.histogram(pcs[:,24],bins)[0]
mhist3 = np.histogram(pcs[:,58],bins)[0]
mhist4 = np.histogram(pcs[:,105],bins)[0]

## PLOTTING

m = np.arange(1000)[::-1]
c = 4e3

v = [0,10,30,100,300,600]

lab = ['#'+str(v[0]+1)+'-'+str(v[1]),\
'#'+str(v[1]+1)+'-'+str(v[2]),\
'#'+str(v[2]+1)+'-'+str(v[3]),\
'#'+str(v[3]+1)+'-'+str(v[4]),\
'#'+str(v[4]+1)+'-'+str(v[5])]

fig,ax = plt.subplots(1)
ax.scatter(s[v[0]:v[1]],k[v[0]:v[1]],30+c*np.sqrt(eigs[v[0]:v[1]]),color=[0.993,  0.906,  0.144],alpha=.8,edgecolors='k',label=lab[0])
ax.scatter(s[v[1]:v[2]],k[v[1]:v[2]],30+c*np.sqrt(eigs[v[1]:v[2]]),color=[ 0.575,  0.844,  0.256],alpha=.8,marker='v',edgecolors='k',label=lab[1])
ax.scatter(s[v[2]:v[3]],k[v[2]:v[3]],30+c*np.sqrt(eigs[v[2]:v[3]]),color=[ 0.246,  0.738,  0.452],alpha=.8,marker='d',edgecolors='k',label=lab[2])
ax.scatter(s[v[3]:v[4]],k[v[3]:v[4]],30+c*np.sqrt(eigs[v[3]:v[4]]),color=[ 0.127,  0.566,  0.550],alpha=.8,marker='p',edgecolors='k',label=lab[3])
ax.scatter(s[v[4]:v[5]],k[v[4]:v[5]],30+c*np.sqrt(eigs[v[4]:v[5]]),color=[ 0.268,  0.009,  0.335],alpha=.8,marker='h',edgecolors='k',label=lab[4])

plt.plot([-2,2],[0,0],'grey')
plt.plot([0,0],[-2,2],'grey')

plt.xlim(-1,1)
plt.ylim(-1,1.5)

plt.legend(loc=4)
plt.xlabel('skewness')
plt.ylabel('kurtosis')
plt.title('Pot. Temp: Principal Components')

### ------------------------

fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
ax1.hist(pcs[:,8],bins,normed=True,alpha=.8,color='g',histtype='step')
ax1.plot(x,gauss,'k')
ax1.set_title('min skewness, mode #9')

ax2.hist(pcs[:,24],bins,normed=True,alpha=.8,color='g',histtype='step')
ax2.plot(x,gauss,'k')
ax2.set_title('max skewness, mode #25')

ax3.hist(pcs[:,58],bins,normed=True,alpha=.8,color='g',histtype='step')
ax3.plot(x,gauss,'k')
ax3.set_title('min kurtosis, mode #59')

ax4.hist(pcs[:,105],bins,normed=True,alpha=.8,color='g',histtype='step')
ax4.plot(x,gauss,'k',label='N(0,1)')
ax4.set_title('max kurtosis, mode #105')
ax4.legend(loc=2)
ax4.set_xlim(-6,6)

## ---------------------------

normal = np.sort(np.random.randn(2778,10000),axis=0)

plt.figure(3)

nq = norm.ppf(np.linspace(1/2778.,1-1/2778.,2778))
nq = norm.ppf(np.linspace(1/2778./2.,1-1/2778./2.,2778))
ra = np.array([-4,4])

plt.plot(nq,np.sort(pcs[:,0],axis=0),'g',lw=.7,label='mode #1-30')
plt.plot(nq,np.sort(pcs[:,range(1,30)],axis=0),'g',lw=.7)
plt.plot(ra,ra,'k')
plt.plot(nq,normal.mean(axis=1),lw=2, label='N(0,1) mean')
plt.fill_between(nq,np.percentile(normal,1,axis=1),np.percentile(normal,99,axis=1),color='grey',alpha=.7)
plt.plot(nq,np.percentile(normal,1,axis=1),'k--')
plt.plot(nq,np.percentile(normal,99,axis=1),'k--',label='99% confidence')
plt.legend(loc=4)
plt.xlabel('Theoretical quantiles')
plt.ylabel('Observed quantiles')
plt.title('Pot. Temp: QQ-Plot of PCs vs N(0,1)')
plt.xlim(-4,4)
plt.ylim(-4,4)
plt.grid(which='both')


plt.show()


