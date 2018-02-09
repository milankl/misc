## WAVE EQ
""" this script solves the 1D-wave equations, i.e.
    
        du/dt = -d/dx(1/2*u^2 + g*h - nu*du/dx) + F
        dh/dt = -d/dx(u*h)."""
    
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time as tictoc
from scipy.signal import tukey

## physical parameters
Lx = 1.     # D = (0,Lx) is the domain
H = 10.      # depth at rest
g = 10.     # gravitational acceleration

## numerical parameters
nx = 200                   # number of grid points for the variable h
nu = nx-1                   # number of grid points for u
nt = 1000                     # number of time steps to integrate

cfl = 0.9                   # cfl number, 1 is the limit

dx = Lx/nx                  # grid spacing

x_h = np.arange(dx/2.,Lx,dx)    # grid x on h-points
x_u = np.arange(dx,Lx-dx+dx/2.,dx)    # on u-points

cph = np.sqrt(g*H)          # long wave phase speed
dt = cfl * dx / cph           # dt based on desired cfl-number

## parameters
rho = 1
v = 0.008
F0 = 500.
F = -F0*np.sin(x_h/Lx*8*np.pi)*np.sin(x_h/Lx*np.pi)**2

## OPERATORS
#GTx = (sparse.diags(np.ones(nx-1),1,shape=(nu,nx)) +\
#   sparse.diags(-np.ones(nx-1),0,shape=(nu,nx))) / dx

GTx = (sparse.diags(np.ones(nx),1,shape=(nu,nx)) +\
    sparse.diags(-np.ones(nx),0,shape=(nu,nx))) / dx
    
Gux = -GTx.T

#LT = Gux.dot(GTx)
#Lu = GTx.dot(Gux)

bc = -0.1
sh = np.zeros(nx)
sh[-1] = bc/dx

def s(t):
    bc = 0.3*np.sin(100*t)
    sh = np.zeros(nx)
    sh[-1] = bc/dx
    return sh

## interpolation from u to T and vice versa
ITu = abs(GTx*dx/2.)
IuT = abs(Gux*dx/2.)

## initial conditions
#h0 = np.zeros(nx)
sigma = 0.08
#h0 = np.exp(-(x_h-.4)**2/(2*sigma**2))
h0 = np.zeros(nx)
u0 = np.zeros(nu)

## time integration - TRIVIAL WAY
h1 = h0.copy()        # store initial conditions
u1 = u0.copy() 

h2 = h0.copy()
u2 = u0.copy()

h3 = h0.copy()
u3 = u0.copy()

h4 = h0.copy()
u4 = u0.copy()


def rhs(u,h):
    h_u = ITu.dot(h+H)
    
    sh = s(t)
    
    # non-linear
    rhs_u = -GTx.dot(.5*(IuT.dot(u**2) + sh*dx/2.) + g*h - v*(Gux.dot(u) + sh))
    rhs_h_dyn = -Gux.dot(u*h_u)
    rhs_h_phy = -sh*(h[-1]+H)  
    
    return rhs_u,rhs_h_dyn,rhs_h_phy

def RK4(u,h):
    k1 = rhs(u,h)
    k2 = rhs(u + dt/2.*k1[0],h + dt/2.*(k1[1]+k1[2]))
    k3 = rhs(u + dt/2.*k2[0],h + dt/2.*(k2[1]+k2[2]))
    k4 = rhs(u + dt*k3[0],h + dt*(k3[1]+k3[2]))
    
    du = (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.
    dh_dyn = (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.
    dh_phy = (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6.
    
    return du,dh_dyn,dh_phy

def SPPT(r):
    #rphase = np.exp(2*np.pi*1j*np.random.rand(nx))
    #r_new = np.real(np.fft.ifft(k0*rphase))    
    #r_new = r_new / r_new.std()
    
    r_new = np.random.randn(1)
    r = c*r + np.sqrt(1-c**2)*r_new
    
    return r

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(9,8))

mass1 = np.empty(nt)
mass2 = np.empty(nt)
mass3 = np.empty(nt)
mass4 = np.empty(nt)

time = np.arange(nt)*dt

a = 0.4
c = 0.5
l = -0.5

k = np.fft.fftfreq(nx,dx)
k0 = np.hstack((0,abs(k[1:])**l)) 

#r = SPPT(np.zeros(nx))
r = SPPT(np.random.randn(1))

t=0

tapering = tukey(2*nx,0.1)[-nx:]

for it in range(nt):
    
    mass1[it] = (h1+H).sum()/(h0+H).sum()
    mass2[it] = (h2+H).sum()/(h0+H).sum()
    mass3[it] = (h3+H).sum()/(h0+H).sum()
    mass4[it] = (h4+H).sum()/(h0+H).sum()
    
    if it == 0:
        l11, = ax1.plot(x_h,h1,'C0',label='deterministic')
        l12, = ax1.plot(x_h,h2,'C1',label='SPPT, BL-tapering, no surf. pert.')
        l13, = ax1.plot(x_h,h3,'C2',label='SPPT, no BL-tapering, no surf. pert.')
        l14, = ax1.plot(x_h,h4,'C3',label='SPPT, no BL-tapering, with surf. pert.')
        
        ax1.set_ylim(8,12)
        ax1.set_xlim(0,1)
        ax1.set_xlabel('x')
        ax1.set_ylabel(r'$h$')
        ax1.set_title(r'surface displacement $h$')
        ax1.legend(loc=2,fontsize=8)
        
        l20, = ax2.plot(time[[0,-1]],[0,0],color='grey',lw=0.5)
        l21, = ax2.plot(time[:it],mass2[:it]-mass1[:it],'C1')
        l22, = ax2.plot(time[:it],mass3[:it]-mass1[:it],'C2')
        l23, = ax2.plot(time[:it],mass4[:it]-mass1[:it],'C3')
        
        ax2.set_ylim(-0.05,.05)
        ax2.set_xlim(0,time[-1])
        ax2.plot([0,nt*dt],[mass1[0],mass1[0]],'grey',lw=0.5)
        ax2.set_xlabel('time')
        ax2.set_ylabel('mass change')
        
        plt.tight_layout()
        
    else:
        r = SPPT(r)
        
        du1,dh_dyn1,dh_phy1 = RK4(u1,h1)
        u1 += dt*du1
        h1 += dt*(dh_dyn1 + dh_phy1)
        
        du2,dh_dyn2,dh_phy2 = RK4(u2,h2)
        u2 += dt*du2
        h2 += dt*dh_dyn2*(1+tapering*a*r.clip(-2,2))
        h2 += dt*dh_phy2
        
        du3,dh_dyn3,dh_phy3 = RK4(u3,h3)
        u3 += dt*du3
        h3 += dt*dh_dyn3*(1+a*r.clip(-2,2))
        h3 += dt*dh_phy3
        
        du4,dh_dyn4,dh_phy4 = RK4(u4,h4)
        u4 += dt*du4
        h4 += dt*dh_dyn4*(1+a*r.clip(-2,2))
        h4 += dt*dh_phy4*(1+a*r.clip(-2,2))
        
        
        t += dt
        
        if it % 2 == 1:
            plt.pause(0.00001)
            l11.set_data(x_h,h1+H)
            
            l12.set_data(x_h,h2+H)            
            l21.set_data(time[:it],mass2[:it]-mass1[:it])
            
            l13.set_data(x_h,h3+H)
            l22.set_data(time[:it],mass3[:it]-mass1[:it])
            
            l14.set_data(x_h,h4+H)
            l23.set_data(time[:it],mass4[:it]-mass1[:it])

plt.pause(2)
plt.close(fig)
#plt.show()        
        


