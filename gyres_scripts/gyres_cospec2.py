## VARIANCE OF HIGH VS LOW

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormaps.py').read())
exec(open('python/ecco2/local_functions.py').read())
from matplotlib.colors import LogNorm
from scipy.signal import detrend

## load data

thi = np.load('python/gyres/temp_highres_sfc.npy')
tlo = np.load('python/gyres/temp_lowres_sfc.npy')
randfield = np.load('python/gyres/patterns/randpattern_ar5_0-100_high.npy')
randfield2 = np.load('python/gyres/patterns/randpattern_ar5_100-200_high.npy')
(time,lat,lon) = np.load('python/gyres/temp_lowres_dim.npy')

## LOAD
(eofL1,pctau1) = np.load('python/gyres/eof_lowres_ltscales.npy')
(eofL2,pctau2) = np.load('python/gyres/eof_highres_ltscales.npy')

(eofs1,pcs1,eigs1) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs2,pcs2,eigs2) = np.load('python/gyres/theta_eofs_highres.npy')

tau1 = np.load('python/gyres/theta_lowres_acftau.npy')
tau2 = np.load('python/gyres/theta_highres_acftau.npy')
print('Data read.')

thi = detrend(thi,axis=0)
tlo = detrend(tlo,axis=0)
print('Data detrended.')

dt = 1.
dy = (lat[1]-lat[0])*111.194
dx = (lon[1]-lon[0])*111.194*np.cos(2*np.pi*lat.mean()/360.)

k,f,thihat = trispec(thi,dt,dy,dx)
print('FFT(T_high) done.')
k,f,tlohat = trispec(tlo,dt,dy,dx)
print('FFT(T_low) done.')
k,f,randhat = trispec(randfield,dt,dy,dx)
print('FFT(T_rand) done.')
k,f,rand2hat = trispec(randfield2,dt,dy,dx)
print('FFT(T_rand) done.')

k = k[1:]
f = f[1:]
thihat = thihat[1:,1:]
tlohat = tlohat[1:,1:]
randhat = randhat[1:,1:]
rand2hat = rand2hat[1:,1:]

v1 = 10.**np.arange(10)

fig1,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
im1 = ax1.contourf(1/k,1/f,tlohat.T,v1,norm = LogNorm())
ax1.contour(1/k,1/f,tlohat.T,v1,colors='k',norm = LogNorm())
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim((1/k).min(),4e3)
ax1.set_ylim((1/f).min(),3e3)
ax1.set_title('T_low')
ax1.set_ylabel('time scale [days]')
ax1.set_xlabel('length scale [km]')

im2 = ax2.contourf(1/k,1/f,thihat.T,v1,norm = LogNorm())
ax2.contour(1/k,1/f,thihat.T,v1,colors='k',norm = LogNorm())

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim((1/k).min(),4e3)
ax2.set_ylim((1/f).min(),3e3)
ax2.set_title('T_high')
ax2.set_xlabel('length scale [km]')

im3 = ax3.contourf(1/k,1/f,randhat.T,v1,norm = LogNorm())
ax3.contour(1/k,1/f,randhat.T,v1,colors='k',norm = LogNorm())
plt.colorbar(im3,ax=(ax1,ax2,ax3))
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_xlim((1/k).min(),4e3)
ax3.set_ylim((1/f).min(),3e3)
ax3.set_title('T_rand')
ax3.set_xlabel('length  scale [km]')

v2 = 10.**np.arange(-2,7)

im4 = ax4.contourf(1/k,1/f,thihat.T/tlohat.T,v2,norm = LogNorm(),cmap=magma)
ax4.contour(1/k,1/f,thihat.T/tlohat.T,v2,colors='k',norm = LogNorm())
ax4.set_yscale('log')
ax4.set_xscale('log')
plt.colorbar(im4,ax=ax4)
ax4.set_xlim((1/k).min(),4e3)
ax4.set_ylim((1/f).min(),3e3)
ax4.set_title('T_high/T_low')
ax4.set_xlabel('length  scale [km]')

## ADD EOF MODES

modemax = 600
v = [0,10,30,100,300,modemax]
s = 4e3
a = 111.194

lab1 = ['#'+str(v[0]+1)+'-'+str(v[1])+'       '+str(round(eigs1[v[0]:v[1]].sum()*100,1))+'% Var',\
'#'+str(v[1]+1)+'-'+str(v[2])+'     '+str(round(eigs1[v[1]:v[2]].sum()*100,1))+'% Var',\
'#'+str(v[2]+1)+'-'+str(v[3])+'   '+str(round(eigs1[v[2]:v[3]].sum()*100,1))+'% Var',\
'#'+str(v[3]+1)+'-'+str(v[4])+'   '+str(round(eigs1[v[3]:v[4]].sum()*100,1))+'% Var',\
'#'+str(v[4]+1)+'-'+str(v[5])+'   '+str(round(eigs1[v[4]:v[5]].sum()*100,1))+'% Var']

lab2 = ['#'+str(v[0]+1)+'-'+str(v[1])+'       '+str(round(eigs2[v[0]:v[1]].sum()*100,1))+'% Var',\
'#'+str(v[1]+1)+'-'+str(v[2])+'     '+str(round(eigs2[v[1]:v[2]].sum()*100,1))+'% Var',\
'#'+str(v[2]+1)+'-'+str(v[3])+'   '+str(round(eigs2[v[2]:v[3]].sum()*100,1))+'% Var',\
'#'+str(v[3]+1)+'-'+str(v[4])+' '+str(round(eigs2[v[3]:v[4]].sum()*100,1))+'% Var',\
'#'+str(v[4]+1)+'-'+str(v[5])+'   '+str(round(eigs2[v[4]:v[5]].sum()*100,1))+'% Var']

alp = .8
zord = 3

#fig,(ax1,ax2) = plt.subplots(2,sharex=True,figsize=(14,10))
ax1.scatter(eofL1[v[0]:v[1]],pctau1[v[0]:v[1]],30+s*np.sqrt(eigs1[v[0]:v[1]]),color=[0.993,  0.906,  0.144],alpha=alp,edgecolors='k',label=lab1[0],zorder=zord)
ax1.scatter(eofL1[v[1]:v[2]],pctau1[v[1]:v[2]],30+s*np.sqrt(eigs1[v[1]:v[2]]),color=[ 0.575,  0.844,  0.256],alpha=alp,marker='v',edgecolors='k',label=lab1[1],zorder=zord)
ax1.scatter(eofL1[v[2]:v[3]],pctau1[v[2]:v[3]],30+s*np.sqrt(eigs1[v[2]:v[3]]),color=[ 0.246,  0.738,  0.452],alpha=alp,marker='d',edgecolors='k',label=lab1[2],zorder=zord)
ax1.scatter(eofL1[v[3]:v[4]],pctau1[v[3]:v[4]],30+s*np.sqrt(eigs1[v[3]:v[4]]),color=[ 0.127,  0.566,  0.550],alpha=alp,marker='p',edgecolors='k',label=lab1[3],zorder=zord)
ax1.scatter(eofL1[v[4]:v[5]],pctau1[v[4]:v[5]],30+s*np.sqrt(eigs1[v[4]:v[5]]),color=[ 0.268,  0.009,  0.335],alpha=alp,marker='h',edgecolors='k',label=lab1[4],zorder=zord)

ax2.scatter(eofL2[v[0]:v[1]],pctau2[v[0]:v[1]],30+s*np.sqrt(eigs2[v[0]:v[1]]),color=[0.993,  0.906,  0.144],alpha=alp,edgecolors='k',label=lab2[0],zorder=zord)
ax2.scatter(eofL2[v[1]:v[2]],pctau2[v[1]:v[2]],30+s*np.sqrt(eigs2[v[1]:v[2]]),color=[ 0.575,  0.844,  0.256],alpha=alp,marker='v',edgecolors='k',label=lab2[1],zorder=zord)
ax2.scatter(eofL2[v[2]:v[3]],pctau2[v[2]:v[3]],30+s*np.sqrt(eigs2[v[2]:v[3]]),color=[ 0.246,  0.738,  0.452],alpha=alp,marker='d',edgecolors='k',label=lab2[2],zorder=zord)
ax2.scatter(eofL2[v[3]:v[4]],pctau2[v[3]:v[4]],30+s*np.sqrt(eigs2[v[3]:v[4]]),color=[ 0.127,  0.566,  0.550],alpha=alp,marker='p',edgecolors='k',label=lab2[3],zorder=zord)
ax2.scatter(eofL2[v[4]:v[5]],pctau2[v[4]:v[5]],30+s*np.sqrt(eigs2[v[4]:v[5]]),color=[ 0.268,  0.009,  0.335],alpha=alp,marker='h',edgecolors='k',label=lab2[4],zorder=zord)

ax1.scatter(eofL1[:modemax],pctau1[:modemax],1,color='k',zorder=zord+1)
ax2.scatter(eofL2[:modemax],pctau2[:modemax],1,color='k',zorder=zord+1)
ax1.legend(loc=1)
ax2.legend(loc=1)

plt.show()



