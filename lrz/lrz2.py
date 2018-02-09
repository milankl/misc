import numpy as np
import matplotlib.pyplot as plt

def lz(x, y, z, s=10, r=28, b=8./3.) :
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot


dt = 0.005
stepCnt = 130

tt = np.linspace(1,stepCnt,stepCnt)

n = 30

# Need one more for the initial values
xs = np.empty((stepCnt + 1,n))
ys = np.empty((stepCnt + 1,n))
zs = np.empty((stepCnt + 1,n))

#Xt = np.empty((stepCnt + 1,3))


# Setting initial values
for i in range(n-2):
    e = np.random.randn(3) * 1e-2
    #Ini3
    #xs[0,i], ys[0,i], zs[0,i] = (8.+e[0] , 8.53+e[1], 28.+e[2])
    xs[0,i], ys[0,i], zs[0,i] = (0.005+e[0] , 0+e[1], 20.+e[2])

ee = np.random.randn(3) * .8e-2
xs[0,-2], ys[0,-2], zs[0,-2] = (0.005+ee[0] , 0+ee[1], 20.+ee[2])
xs[0,-1], ys[0,-1], zs[0,-1] = (0.005, 0., 20.)
    
#and for the "truth"
#Xt[0,0], Xt[0,1], Xt[0,2] = (8,8.53,28)
#Xt[0,0], Xt[0,1], Xt[0,2] = (8,8.53,28)

# Stepping through "time".
for i in range(stepCnt) :
    # Derivatives of the X, Y, Z state
    
    #euler forward
    dxyz = np.empty((3,n-1))
    txyz = np.empty(3)
    txyz2 = np.empty(3)
    txs = np.empty(1)
    tys = np.empty(1)
    tzs = np.empty(1)
    
    for ii in range(n-1):
        dxyz[0,ii], dxyz[1,ii], dxyz[2,ii] = lz(xs[i,ii], ys[i,ii], zs[i,ii])
    
    xs[i+1,:-1] = xs[i,:-1] + (dxyz[0,:] * dt)
    ys[i+1,:-1] = ys[i,:-1] + (dxyz[1,:] * dt)
    zs[i+1,:-1] = zs[i,:-1] + (dxyz[2,:] * dt)
    
    txyz[0], txyz[1], txyz[2] = lz(xs[i,-1], ys[i,-1], zs[i,-1])
    
    txs = xs[i,-1] + (txyz[0] * dt * (1-1e-1))
    tys = ys[i,-1] + (txyz[1] * dt * (1-1e-1))
    tzs = zs[i,-1] + (txyz[2] * dt * (1-1e-1))
    
    txyz2[0], txyz2[1], txyz2[2] = lz(txs, tys, tzs)
    
    xs[i+1,-1] = txs + (txyz2[0] * dt * 1e-1)
    ys[i+1,-1] = tys + (txyz2[1] * dt * 1e-1)
    zs[i+1,-1] = tzs + (txyz2[2] * dt * 1e-1)    
    
    
#     heun for the truth
#     #txyz1 = np.empty(3)
#     txyz2 = np.empty(3)
#     txyz1[0], txyz1[1], txyz1[2] = lz(Xt[i,0],Xt[i,1],Xt[i,2])
#     Xt1 = Xt[i,:] + txyz1 * dt
#     #corrector
#     txyz2[0], txyz2[1], txyz2[2] = lz(Xt1[0],Xt1[1],Xt1[2])
#     Xt[i+1,:] = Xt[i,:] + .5 * dt * (txyz1 + txyz2)

mxs = np.mean(xs,1)

#plotting
fig, ax = plt.subplots()

line = []
poi = []

for t in range(stepCnt):
    if t == 0:
        for ii in range(n-2):
            line1,poi1 = ax.plot(tt[0],xs[0,ii],'k',tt[0],xs[0,ii],'k.',lw=0.6,ms=13)
            line.append(line1)
            poi.append(poi1) 
        
        #best estimate
        line1,poi1 = ax.plot(tt[0],xs[0,-2],'b-',tt[0],xs[0,-2],'b.',lw=1,ms=20)
        line.append(line1)
        poi.append(poi1) 
        #truth
        line1,poi1 = ax.plot(tt[0],xs[0,-1],'g-',tt[0],xs[0,-2],'g.',lw=1,ms=20)
        line.append(line1)
        poi.append(poi1) 
        #mean
        line1,poi1 = ax.plot(tt[0],mxs[0],'r-',tt[0],mxs[0],'r.',lw=1,ms=20)
        line.append(line1)
        poi.append(poi1) 
        
        ax.set_xlim(10, 140) 
        ax.set_ylim(-3, 3)
        ax.set_xlabel('time [steps]')
        ax.set_ylabel('x')
        ax.text(30,2,'truth',bbox=dict(facecolor='green', alpha=0.5),fontsize=15)
        ax.text(30,1.7,'best estimate',bbox=dict(facecolor='blue', alpha=0.5),fontsize=15)
        ax.text(30,1.4,'ensemble member',bbox=dict(facecolor='black', alpha=0.5),fontsize=15)
        ax.text(30,1.1,'ensemble mean',bbox=dict(facecolor='red', alpha=0.5),fontsize=15)
        
    else:
        ii = 0
        for l in line:
            if ii == n:
                l.set_data(tt[:t], mxs[:t])
            else:
                l.set_data(tt[:t], xs[:t,ii])
            ii = ii+1
        
        ii = 0
        for p in poi:
            if ii == n:
                p.set_data(tt[t], mxs[t])
            else:
                p.set_data(tt[t], xs[t,ii])
            ii = ii+1
        
        
        ax.set_title('t='+str(t)+' / '+str(stepCnt))
        
    plt.pause(0.001)
    t=t+1
plt.pause(15)
plt.close('all')
