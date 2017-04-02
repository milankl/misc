## ANALYSE EOF TIME AND SPACE SCALES

import numpy as np
import matplotlib.pyplot as plt
exec(open('python/ecco2/colormap.py').read())

## LOAD
(eofL,pctau) = np.load('python/ecco2/eof_ltscales.npy')
(eofs,pcs,eigs) = np.load('python/ecco2/theta_eofs_1000.npy')

acftau = np.load('python/ecco2/theta_acftau_detr.npy')
lscale,m = np.load('python/ecco2/decorr_lscale_1.npy')

eccomask = np.load('python/ecco2/ecco_mask_NA.npy')

## apply mask
a = 111.194
perc = 2

acftau = np.ma.masked_array(acftau,mask=eccomask)
lscale = np.ma.masked_array(lscale,m)

eccoL = [np.percentile(lscale[lscale > 0].data,perc),lscale.max()]
eccotau = [acftau.min(),acftau.max()]

modemax = 600
modes = np.arange(modemax)

eofL = eofL/2

modei = 0
modej = 100

x = eofL[modei:modej].min()
y = pctau[modei:modej].min()
w = eofL[modei:modej].max() - x
h = pctau[modei:modej].max() - y

##

xe,ye = eccoL[0],eccotau[0]
we,he = eccoL[1]-xe,eccotau[1]-ye

##

modei2 = 10
modej2 = 40

x2 = eofL[modei2:modej2].min()
y2 = pctau[modei2:modej2].min()
w2 = eofL[modei2:modej2].max() - x2
h2 = pctau[modei2:modej2].max() - y2

##

modei3 = 30
modej3 = 60

x3 = eofL[modei3:modej3].min()
y3 = pctau[modei3:modej3].min()
w3 = eofL[modei3:modej3].max() - x3
h3 = pctau[modei3:modej3].max() - y3


## GET SCALES FROM DATA

path = 'python/ecco2/patterns/'

t1 = np.ma.masked_array(np.load(path+'randpattern_1_acftau.npy'),mask=eccomask)
t2 = np.ma.masked_array(np.load(path+'randpattern_2_acftau.npy'),mask=eccomask)
t3 = np.ma.masked_array(np.load(path+'randpattern_3_acftau.npy'),mask=eccomask)

l1 = np.ma.masked_array(np.load(path+'decorr_lscale_1.npy')[0,:,:],mask=m)
l2 = np.ma.masked_array(np.load(path+'decorr_lscale_2.npy')[0,:,:],mask=m)
l3 = np.ma.masked_array(np.load(path+'decorr_lscale_3.npy')[0,:,:],mask=m)

l1min = np.percentile(l1[l1 > 0].data,perc)
l2min = np.percentile(l2[l2 > 0].data,perc)
l3min = np.percentile(l3[l3 > 0].data,perc)

## PLOTTING

v = [0,10,30,100,300,modemax]
s = 4e3

lab = ['#'+str(v[0]+1)+'-'+str(v[1])+'       '+str(round(eigs[v[0]:v[1]].sum()*100,1))+'% Var',\
'#'+str(v[1]+1)+'-'+str(v[2])+'     '+str(round(eigs[v[1]:v[2]].sum()*100,1))+'% Var',\
'#'+str(v[2]+1)+'-'+str(v[3])+'   '+str(round(eigs[v[2]:v[3]].sum()*100,1))+'% Var',\
'#'+str(v[3]+1)+'-'+str(v[4])+' '+str(round(eigs[v[3]:v[4]].sum()*100,1))+'% Var',\
'#'+str(v[4]+1)+'-'+str(v[5])+'   '+str(round(eigs[v[4]:v[5]].sum()*100,1))+'% Var']

fig,ax = plt.subplots(1)
ax.scatter(eofL[v[0]:v[1]],pctau[v[0]:v[1]],30+s*np.sqrt(eigs[v[0]:v[1]]),color=[0.993,  0.906,  0.144],alpha=.8,edgecolors='k',label=lab[0])
ax.scatter(eofL[v[1]:v[2]],pctau[v[1]:v[2]],30+s*np.sqrt(eigs[v[1]:v[2]]),color=[ 0.575,  0.844,  0.256],alpha=.8,marker='v',edgecolors='k',label=lab[1])
ax.scatter(eofL[v[2]:v[3]],pctau[v[2]:v[3]],30+s*np.sqrt(eigs[v[2]:v[3]]),color=[ 0.246,  0.738,  0.452],alpha=.8,marker='d',edgecolors='k',label=lab[2])
ax.scatter(eofL[v[3]:v[4]],pctau[v[3]:v[4]],30+s*np.sqrt(eigs[v[3]:v[4]]),color=[ 0.127,  0.566,  0.550],alpha=.8,marker='p',edgecolors='k',label=lab[3])
ax.scatter(eofL[v[4]:v[5]],pctau[v[4]:v[5]],30+s*np.sqrt(eigs[v[4]:v[5]]),color=[ 0.268,  0.009,  0.335],alpha=.8,marker='h',edgecolors='k',label=lab[4])

ax.scatter(eofL[:modemax],pctau[:modemax],1,color='k')
ax.legend(loc=4)

ax.add_patch(plt.Rectangle([x,y],w,h,fill=False,edgecolor='k',alpha=.5))
ax.add_patch(plt.Rectangle([x2,y2],w2,h2,fill=False,edgecolor='k',alpha=.5))
ax.add_patch(plt.Rectangle([x3,y3],w3,h3,fill=False,edgecolor='k',alpha=.5))
ax.add_patch(plt.Rectangle([xe,ye],we,he,color='g',alpha=.1))

ax.text(xe,ye,'ECCO2',verticalalignment='bottom')
ax.text(x,y+h,'MODE '+str(modei+1)+'-'+str(modej),verticalalignment='bottom')
ax.text(x2,y2+h2,'MODE '+str(modei2+1)+'-'+str(modej2),verticalalignment='bottom')
ax.text(x3,y3+h3,'MODE '+str(modei3+1)+'-'+str(modej3),verticalalignment='bottom')

#reconstruct 1-100
ax.add_patch(plt.Rectangle([l1min,t1.min()],l1.max()-l1min,t1.max()-t1.min(),alpha=.15,facecolor='y'))
ax.add_patch(plt.Rectangle([l2min,t2.min()],l2.max()-l2min,t2.max()-t2.min(),alpha=.15,facecolor='g'))
ax.add_patch(plt.Rectangle([l3min,t3.min()],l3.max()-l3min,t3.max()-t3.min(),alpha=.15,facecolor='c'))

ax.text(l1.max(),t1.min(),'RANDOM PATTERN '+str(modei+1)+'-'+str(modej),verticalalignment='bottom',\
horizontalalignment='right')
ax.text(l2.max(),t2.min(),'RANDOM PATTERN '+str(modei2+1)+'-'+str(modej2),verticalalignment='bottom',\
horizontalalignment='right')
ax.text(l3.max(),t3.min(),'RANDOM PATTERN '+str(modei3+1)+'-'+str(modej3),verticalalignment='bottom',\
horizontalalignment='right')

ax.plot([a,a],[1,1000],'b--')
ax.plot([a/4,a/4],[1,1000],'g--')
ax.plot([1,10e4],[3,3],'g--')
ax.text(a+1,1.2,r'1$\degree{}$',rotation='vertical')
ax.text(a/4+1,1.2,r'1/4$\degree{}$',rotation='vertical')
ax.text(21,3.1,'3 days')

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(20,5000)
ax.set_ylim(1,500)
ax.set_xlabel('Length scale [km]')
ax.set_ylabel('Time scale [days]')
ax.set_title('Pot. Temp: Length and time scales')
plt.show()
