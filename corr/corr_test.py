## Correlation increase with running mean
import numpy as np
import matplotlib.pyplot as plt
from cmocean import cm

nt = 150    # length of time series
ni = 100   # number of samples

arclist = np.arange(0.05,1,.1)
trueclist = np.arange(0,1,.1)

narc = len(arclist)
ntruec = len(trueclist)

corrm = np.empty(ni)
corrd = np.empty((narc,ntruec))
trueclistf = np.empty_like(corrd)
arclistf = np.empty_like(corrd)

def rmean(x, N):
    """ cutting off the edges. """
    s = int(N-1)
    return np.convolve(x, np.ones((N,))/N)[s:-s]

def ar1(n,arc):
    x = np.empty(n)
    x[0] = np.random.randn(1)
    for i in range(n-1):
        x[i+1] = arc*x[i] + np.sqrt(1-arc**2)*np.random.randn(1)
    
    return x

def acf(x,l):
    """ autocorrelation function of vector x up to lag l."""
    return np.array([1]+[np.corrcoef(x[:-i],x[i:])[0,1] for i in range(1,l)])
    
for l in range(ntruec):
    c = trueclist[l]
    print('True correlation is %.2f' % c)
        
    for k in range(narc):
        arc = arclist[k]
        print('Autocorrelation is %.2f' % arc)
        
        rl = 2*int(np.ceil(-1/np.log(arc)))
        nt = 100*rl
        print(rl)
    
        for i in range(ni):
            x = ar1(nt,arc)
            y = c*x + np.sqrt(1-c**2)*np.random.randn(nt)
            
            xm = x.reshape((100,rl)).mean(axis=1)
            ym = y.reshape((100,rl)).mean(axis=1)
            
            corrm[i] = np.corrcoef(xm,ym)[0,1]
        
        corrd[k,l] = (corrm-c).mean()
        trueclistf[k,l] = c
        arclistf[k,l] = arc
            
## plotting
fig,ax = plt.subplots(1,1)

corrdf = (corrd+trueclistf).flatten()
trueclistfs = trueclistf.flatten()
arclistfs = arclistf.flatten()

thermal_discrete = cm.thermal.from_list('thermal_discrete',cm.thermal(np.linspace(0,1,10)),10)

q = ax.scatter(trueclistfs,corrdf,50,arclistfs,cmap=thermal_discrete)
cb = plt.colorbar(q,ticks=arclist,drawedges=True)
cb.set_label('AR1-coefficient')
q.set_clim(0.,1.)
ax.plot([-.2,1],[-.2,1],'k--',lw=2)
ax.set_xlim(-.1,1)
ax.set_ylim(-.1,1)
ax.set_xlabel('True correlation')
ax.set_ylabel('Correlation after rmean')
ax.set_title('Correlation increase for running mean filtered time series')
plt.tight_layout()
plt.show()
