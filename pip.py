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


# Define the polygon
poly_coords = [[0.3,0.3],[.1,0.4],[.4,.6],[.7,.5],[.6,.2]]
poly_path = path.Path(poly_coords)

# create points and check whether they are in polygon
x = np.linspace(0,1,50)
xx,yy = np.meshgrid(x,x)
flags = poly_path.contains_points(random_points)

# plotting
fig,ax = plt.subplots(1,1)

p = PatchCollection([Polygon(poly_coords)], alpha=0.4)
ax.add_collection(p)

ax.scatter(random_points[:,0],random_points[:,1])


plt.show() 