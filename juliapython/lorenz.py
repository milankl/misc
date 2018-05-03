from time import time as tictoc

def time_integration(N,x,y,z):
    
    # pre-allocate
    dx,dy,dz = 0.,0.,0.
    
    for i in range(N-1):
        
        # RHS
        dx = dt*(sig*(y-x))
        dy = dt*(x*(rho-z) - y)
        dz = dt*(x*y - beta*z)
        
        # Euler forward
        x += dt*dx
        y += dt*dy
        z += dt*dz
    
    return x,y,z


N = int(1e8)        # number of time steps
dt = 0.01   # time step

# initial conditions
x,y,z = 5.,5.,5.    

# L63 parameters
sig = 10.
rho = 28.
beta = 8./3.


# RUN
t0 = tictoc()
x,y,z = time_integration(N,x,y,z)
t1 = tictoc()-t0
print("time elapsed {:.2f}s".format(t1))

print((x,y,z))

