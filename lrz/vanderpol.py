# illustrate the van der pol oscillator

import numpy as np
import matplotlib.pyplot as plt

def rhs(x, y, mu=1.):
    xp = mu*(y - (x**3)/3. + x)
    yp = -1./mu * x
    return xp, yp


dt = 1e-2
steps = 1000
t = np.arange(0,steps*dt,dt)

s = 20
n = s**2
l = 3

# Need one more for the initial values
xs = np.nan*np.empty((steps + 1,n))
ys = np.nan*np.empty((steps + 1,n))

x = np.linspace(-l,l,s)
y = np.linspace(-l,l,s)

xx,yy = np.meshgrid(x,y)

# Setting initial values
xs[0,:] = xx.reshape(n)
ys[0,:] = yy.reshape(n)


# Stepping
for i in range(steps):
    # Derivatives of the X, Y, Z state
    
    #heun
    dxy = np.empty((2,n))
    dxyc = np.empty((2,n)) #corrector
    
    dxy[0,:], dxy[1,:] = rhs(xs[i,:], ys[i,:])
    
    xs[i+1,:] = xs[i,:] + dt*dxy[0,:]
    ys[i+1,:] = ys[i,:] + dt*dxy[1,:]
    
    dxyc[0,:], dxyc[1,:] = rhs(xs[i+1,:], ys[i+1,:])
    
    xs[i+1,:] = xs[i,:] + .5*dt*(dxy[0,:] + dxyc[0,:])
    ys[i+1,:] = ys[i,:] + .5*dt*(dxy[1,:] + dxyc[1,:])    

##plotting
plt.figure(1)
plt.plot(xs[:300,:],ys[:300,:],'black',xs[-300:,:],ys[-300:,:],'blue',lw=.5)
plt.xlim(-l,l)
plt.ylim(-l,l)
plt.xlabel('x')
plt.ylabel('y')

vx = (xs[1,:]-xs[0,:]).reshape(s,s)
vy = (ys[1,:]-ys[0,:]).reshape(s,s)

plt.figure(2)
plt.quiver(xx,yy,vx,vy)
plt.xlim(-l,l)
plt.ylim(-l,l)
plt.xlabel('x')
plt.ylabel('y')
plt.show()