## ANALYSE EOF TIME AND SPACE SCALES

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormap.py').read())

## LOAD
(eofL1,pctau1) = np.load('python/gyres/eof_lowres_ltscales.npy')
(eofL2,pctau2) = np.load('python/gyres/eof_highres_ltscales.npy')

(eofs1,pcs1,eigs1) = np.load('python/gyres/theta_eofs_lowres.npy')
(eofs2,pcs2,eigs2) = np.load('python/gyres/theta_eofs_highres.npy')

tau1 = np.load('python/gyres/theta_lowres_acftau.npy')
tau2 = np.load('python/gyres/theta_highres_acftau.npy')

#lscale,m = np.load('python/ecco2/decorr_lscale_1.npy')

## PLOTTING

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

fig,(ax1,ax2) = plt.subplots(2,sharex=True,figsize=(14,10))
ax1.scatter(eofL1[v[0]:v[1]],pctau1[v[0]:v[1]],30+s*np.sqrt(eigs1[v[0]:v[1]]),color=[0.993,  0.906,  0.144],alpha=.8,edgecolors='k',label=lab1[0])
ax1.scatter(eofL1[v[1]:v[2]],pctau1[v[1]:v[2]],30+s*np.sqrt(eigs1[v[1]:v[2]]),color=[ 0.575,  0.844,  0.256],alpha=.8,marker='v',edgecolors='k',label=lab1[1])
ax1.scatter(eofL1[v[2]:v[3]],pctau1[v[2]:v[3]],30+s*np.sqrt(eigs1[v[2]:v[3]]),color=[ 0.246,  0.738,  0.452],alpha=.8,marker='d',edgecolors='k',label=lab1[2])
ax1.scatter(eofL1[v[3]:v[4]],pctau1[v[3]:v[4]],30+s*np.sqrt(eigs1[v[3]:v[4]]),color=[ 0.127,  0.566,  0.550],alpha=.8,marker='p',edgecolors='k',label=lab1[3])
ax1.scatter(eofL1[v[4]:v[5]],pctau1[v[4]:v[5]],30+s*np.sqrt(eigs1[v[4]:v[5]]),color=[ 0.268,  0.009,  0.335],alpha=.8,marker='h',edgecolors='k',label=lab1[4])

ax2.scatter(eofL2[v[0]:v[1]],pctau2[v[0]:v[1]],30+s*np.sqrt(eigs2[v[0]:v[1]]),color=[0.993,  0.906,  0.144],alpha=.8,edgecolors='k',label=lab2[0])
ax2.scatter(eofL2[v[1]:v[2]],pctau2[v[1]:v[2]],30+s*np.sqrt(eigs2[v[1]:v[2]]),color=[ 0.575,  0.844,  0.256],alpha=.8,marker='v',edgecolors='k',label=lab2[1])
ax2.scatter(eofL2[v[2]:v[3]],pctau2[v[2]:v[3]],30+s*np.sqrt(eigs2[v[2]:v[3]]),color=[ 0.246,  0.738,  0.452],alpha=.8,marker='d',edgecolors='k',label=lab2[2])
ax2.scatter(eofL2[v[3]:v[4]],pctau2[v[3]:v[4]],30+s*np.sqrt(eigs2[v[3]:v[4]]),color=[ 0.127,  0.566,  0.550],alpha=.8,marker='p',edgecolors='k',label=lab2[3])
ax2.scatter(eofL2[v[4]:v[5]],pctau2[v[4]:v[5]],30+s*np.sqrt(eigs2[v[4]:v[5]]),color=[ 0.268,  0.009,  0.335],alpha=.8,marker='h',edgecolors='k',label=lab2[4])

ax1.scatter(eofL1[:modemax],pctau1[:modemax],1,color='k')
ax2.scatter(eofL2[:modemax],pctau2[:modemax],1,color='k')
ax1.legend(loc=4)
ax2.legend(loc=1)

ax1.plot([a,a],[1,1000],'b--')
ax1.plot([a/4,a/4],[1,1000],'g--')
ax1.plot([1,10e4],[1,1],'g--')
ax1.text(a-8,80,r'1$\degree{}$',rotation='vertical')
ax1.text(a/4-2,80,r'1/4$\degree{}$',rotation='vertical')
ax1.text(21,1.1,'1 day')

ax2.plot([a,a],[1,1000],'b--')
ax2.plot([a/4,a/4],[1,1000],'g--')
ax2.plot([1,10e4],[1,1],'g--')
ax2.text(a-8,80,r'1$\degree{}$',rotation='vertical')
ax2.text(a/4-2,80,r'1/4$\degree{}$',rotation='vertical')
ax2.text(21,1.1,'1 day')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim(20,1000)
ax1.set_ylim(.5,100)
ax1.set_ylabel('Time scale [days]')
ax1.set_title('T_low, length and time scales')

ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xlim(20,1000)
ax2.set_ylim(.5,100)
ax2.set_xlabel('Length scale [km]')
ax2.set_ylabel('Time scale [days]')
ax2.set_title('T_high, length and time scales')
plt.show()
