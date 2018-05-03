import numpy as np
from matplotlib import path
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

"""
first = -3
size  = (3-first)/100
xv,yv = np.meshgrid(np.linspace(-3,3,100),np.linspace(-3,3,100))
p = path.Path([(0,0), (0, 1), (1, 1), (1, 0)])  # square with legs length 1 and bottom left corner at the origin
flags = p.contains_points(np.hstack((xv.flatten()[:,np.newaxis],yv.flatten()[:,np.newaxis])))


grid = np.zeros((101,101),dtype='bool')
grid[((xv.flatten()-first)/size).astype('int'),((yv.flatten()-first)/size).astype('int')] = flags

xi,yi = np.random.randint(-300,300,100)/100,np.random.randint(-300,300,100)/100
vflag = grid[((xi-first)/size).astype('int'),((yi-first)/size).astype('int')]

#plt.imshow(grid.T,origin='lower',interpolation='nearest',cmap='binary')
plt.scatter(((xi-first)/size).astype('int'),((yi-first)/size).astype('int'),c=vflag,cmap='Greens',s=90)
plt.show()

"""

polygon = Polygon(np.random.rand(N, 2), True)
patches.append(polygon)

colors = 100*np.random.rand(len(patches))
p = PatchCollection(patches, alpha=0.4)
p.set_array(np.array(colors))
ax.add_collection(p)
fig.colorbar(p, ax=ax)

plt.show() 