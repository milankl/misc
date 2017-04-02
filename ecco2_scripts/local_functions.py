## LOCAL FUNCTIONS

import numpy as np

def acf(x,l):
    """ autocorrelation function of vector x up to lag l."""
    return np.array([1]+[np.corrcoef(x[:-i],x[i:])[0,1] for i in range(1,l)])

def acfast(x,l):
    """ autocorrelation function of vector x up to lag l. Assumes zero mean for x,
    which does not change over time. Assumes a constant variance of x over time.
    Only appropriate for l << length(x)."""
    v = np.dot(x,x)/len(x)
    return np.array([1]+[np.dot(x[:-i],x[i:])/(len(x)-i)/v for i in range(1,l)])

def findzc(x,a):
    """ find crossing of vector x with value a, several values might be returned."""
    return np.where(abs(np.diff(np.ceil(x-a))) == 1)[0]

def mcorr(x,y):
    """correlation between vector x and 2d-array y. faster than np.corrcoef. x,y are assumed to have zero mean in time dimension (first dim) and variance = 1."""
    yr = np.ma.reshape(y,(y.shape[0],np.prod(y.shape[1:])))
    return (np.ma.dot(x,yr).reshape(y.shape[1:]) / (x.shape[0] - 1))
    

def trispec(a,dt,dx,dy):
    """ Computes a wavenumber-frequency plot for 3D (t,x,y) data via radial (k = sqrt(kx**2 + ky**2)) integration. TODO: correct normalisation, so that the integral in normal space corresponds to the integral in Fourier space.
    """
    
    nt,ny,nx = np.shape(a)
    kx = (1/(dx))*np.hstack((np.arange(0,(nx+1)/2.),np.arange(-nx/2.+1,0)))/float(nx)
    ky = (1/(dy))*np.hstack((np.arange(0,(ny+1)/2.),np.arange(-ny/2.+1,0)))/float(ny)
    f = (1/(dt))*np.hstack((np.arange(0,(nt+1)/2.),np.arange(-nt/2.+1,0)))/float(nt)

    kxx,kyy = np.meshgrid(kx,ky)
    # radial distance from kx,ky = 0
    kk = np.sqrt(kxx**2 + kyy**2) 

    if nx >= ny: #kill negative wavenumbers
        k  = kx[:int(nx/2)+1]
    else:
        k  = ky[:int(ny/2)+1]

    f = f[:int(nt/2)+1] #kill negative frequencies
    dk = k[1] - k[0]

    # create radial coordinates, associated with k[i]
    # nearest point interpolation to get points within the -.5,.5 annulus
    rcoords = []
    for i in range(len(k)):
        rcoords.append(np.where((kk>(k[i]-.5*dk))*(kk<=(k[i]+.5*dk))))

    # 3D FFT
    p = np.fft.fftn(a)
    p = np.real(p * np.conj(p))

    # mulitply by dk to have the corresponding integral
    spec = np.zeros((len(k),len(f)))
    for i in range(len(k)):
        spec[i,:] = np.sum(p[:int(nt/2)+1,rcoords[i][0],rcoords[i][1]],axis=1)*dk

    return k,f,spec
