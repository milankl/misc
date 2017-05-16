## 2D nearest neighbour search
import numpy as np
import matplotlib.pyplot as plt

x0,x1 = 0,1     # domain edges
y0,y1 = 0,1

N = 50        # number of points in the domain (excluding the origin)

Xi = x0 + x1*np.random.rand(1)  # one random point as origin
Yi = y0 + y1*np.random.rand(1)

X = x0 + x1*np.random.rand(N)   # X,Y positions of all other points
Y = y0 + y1*np.random.rand(N)

# find nearest neighbour to point i
nn = np.argmin((X-Xi)**2+(Y-Yi)**2)

# radius for plotting purposes
r = np.sqrt((Xi-X[nn])**2 + (Yi-Y[nn])**2)

# plotting
fig,ax = plt.subplots(1,figsize=(6,6))

ax.plot(Xi,Yi,'C1o')
ax.plot(X[nn],Y[nn],'C2o')
ax.scatter(X,Y,s=10)

circle = plt.Circle((Xi, Yi), r, alpha=0.5)
ax.add_artist(circle)

ax.set_xlim(x0,x1)
ax.set_ylim(y0,y1)

plt.show()
