import numpy as np
from matplotlib import path
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# Define the polygon
poly_coords = [[0.3,0.3],[.1,0.4],[.4,.6],[.7,.5],[.6,.2]]
poly_path = path.Path(poly_coords)

# create points and check whether they are in polygon
x = np.linspace(0,1,10)
xx,yy = np.meshgrid(x,x)
xx,yy = xx.flatten(),yy.flatten()
flags = poly_path.contains_points(np.vstack((xx,yy)).T)

## PLOTTING
fig,ax = plt.subplots(1,1)

# plot the polygon
p = PatchCollection([Polygon(poly_coords)], alpha=0.4)
ax.add_collection(p)

# plot grid points
ax.scatter(xx[flags],yy[flags],color="C1")
ax.scatter(xx[~flags],yy[~flags],color="C2")

plt.show() 