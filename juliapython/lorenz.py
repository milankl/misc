from time import time as tictoc

def time_integration(N,x,y,z):

    dx,dy,dz = 0.,0.,0              # pre-allocate

    for i in range(N-1):

        dx = dt*(sig*(y-x))         # RHS
        dy = dt*(x*(rho-z) - y)
        dz = dt*(x*y - beta*z)

        x += dt*dx                  # Euler forward
        y += dt*dy
        z += dt*dz

    return x,y,z

N = int(1e8)                        # number of time steps
dt = 0.01                           # time step

x,y,z = 5.,5.,5.                    # initial conditions

sig = 10.                           # L63 parameters
rho = 28.
beta = 8./3.

t0 = tictoc()                       # RUN
x,y,z = time_integration(N,x,y,z)
t1 = tictoc()-t0
print("time elapsed {:.2f}s".format(t1))
