import numpy as np
import matplotlib.pyplot as plt

def lz(x, y, z, s=10, r=28, b=8./3.) :
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot


dt = 0.01
stepCnt = 300

tt = np.linspace(1,stepCnt,stepCnt)

n = 30

# Need one more for the initial values
xs = np.empty((stepCnt + 1,n))
ys = np.empty((stepCnt + 1,n))
zs = np.empty((stepCnt + 1,n))


# Setting initial values
for i in range(n):
    e = np.random.randn(3) * 1e-1
    #Ini3
    xs[0,i], ys[0,i], zs[0,i] = (10.+e[0] , 8.53+e[1], 28.+e[2])
    if i > 15:
        xs[0,i], ys[0,i], zs[0,i] = (-10.+e[0] , 2+e[1], 28.+e[2])
    else:
        xs[0,i], ys[0,i], zs[0,i] = (10.+e[0] , 8.53+e[1], 28.+e[2])
        
# Stepping through "time".
for i in range(stepCnt) :
    # Derivatives of the X, Y, Z state
    
    #euler forward
    dxyz = np.empty((3,n))
    
    for ii in range(n):
        dxyz[0,ii], dxyz[1,ii], dxyz[2,ii] = lz(xs[i,ii], ys[i,ii], zs[i,ii])
    
    xs[i+1,:] = xs[i,:] + (dxyz[0,:] * dt)
    ys[i+1,:] = ys[i,:] + (dxyz[1,:] * dt)
    zs[i+1,:] = zs[i,:] + (dxyz[2,:] * dt)

#plotting
fig, ax = plt.subplots()

line = []
poi = []

for t in range(stepCnt):
    if t == 0:
        for ii in range(n):
            ax.plot(xs[:,ii],zs[:,ii],'grey',lw=0.2)
            line1,poi1 = ax.plot(xs[0,ii],zs[0,ii],'b',xs[0,ii],zs[0,ii],'r.',lw=0.4,ms=10)
            line.append(line1)
            poi.append(poi1) 

        
        ax.set_xlim(-30, 30) 
        ax.set_ylim(0, 55)
        ax.set_xlabel('x')
        ax.set_ylabel('z')
        
    else:
        ii = 0
        for l in line:
            l.set_data(xs[:t,ii], zs[:t,ii])
            ii = ii+1
        
        ii = 0
        for p in poi:
            p.set_data(xs[t,ii], zs[t,ii])
            ii = ii+1
        
        
        ax.set_title('t='+str(t)+' / '+str(stepCnt))
        
    plt.pause(0.0001)
    
    
plt.pause(10)
plt.close('all')
