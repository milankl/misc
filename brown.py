## BROWNIAN MOTION
import numpy as np
import matplotlib.pyplot as plt

## parameters

m = 1.
n = 100
dt = 1e-2
s = 0.03

## initial conditions

N = 2      # number of particles
#r = np.random.rand(N,2)*.9+.05
#uv = (np.random.rand(N,2)-.5)*2e-1+
r = np.array([[0.41,.5],[0.6,.5]]).T
uv = np.array([[.3,0],[-.3,0]]).T


##
global ncol
ncol = 0
#status matrix, True means allowed to collide, False not allowed
X,Y = np.array([r[:,0]]*N),np.array([r[:,1]]*N)
S = np.sqrt((X - X.T)**2 + (Y - Y.T)**2) > s

def border(r,uv):
    x,y = r[:,0],r[:,1]
    u,v = uv[:,0],uv[:,1]
        
    xr = (x > 1)
    xl = (x < 0)
    yu = (y > 1)
    yl = (y < 0)
    
    x[xr] = 2 - x[xr]
    x[xl] = -x[xl]
    u[xr+xl] = -u[xr+xl]
    
    y[yu] = 2 - y[yu]
    y[yl] = -y[yl]
    v[yu+yl] = -v[yu+yl]
    
    return np.vstack((x,y)).T,np.vstack((u,v)).T

def collision(r1,r2,uv1,uv2):
    global ncol
    ncol += 1
    
    dr = r1-r2
    duv = uv1-uv2
    
    uv1 = uv1 - duv.dot(dr) / dr.dot(dr) * dr
    uv2 = uv2 + duv.dot(dr) / dr.dot(dr) * dr
    
    print(ncol,end=',')
    
    return uv1,uv2

for i in range(n-1):

    r = r + uv*dt
    r,uv = border(r,uv)
    
    X,Y = np.array([r[:,0]]*N),np.array([r[:,1]]*N)
    dR = np.sqrt((X - X.T)**2 + (Y - Y.T)**2)
    col_ind = np.argwhere(np.triu((dR < s)*S,1))
    
    for j,k in col_ind:
        uv[j],uv[k] = collision(r[j],r[k],uv[j],uv[k])
    
    # update status matrix
    S = dR > s
    
    if i == 0:
        fig,ax = plt.subplots(1,1)
        q = ax.scatter(r[:,0],r[:,1],50,np.arange(N))
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
    else:
        q.set_offsets(r.T)
        plt.pause(.005)

plt.pause(3)
plt.close(fig)
