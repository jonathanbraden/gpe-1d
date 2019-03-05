#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

#dx = 1.0918300671383692/2.
dx = 0.54591503356928461
dt = 0.17562036827601815
om = 4.4721359549995796
lv = 1.4
nu = 2.e-3

def fourier_derivative(f,dx):
    n = f.size; norm = 2.*np.pi/dx
    fk = 1j*np.fft.rfft(f)*np.fft.rfftfreq(n)*norm
    return np.fft.irfft(fk)

def fourier_laplacian(f,dx):
    return

def fourier_laplacian_2deriv(f,dx):
    return

def fourier_derivative_cmpx(f,dx):
    return

# These functions assume that the variables are being stored in a single array with the last column indexing the fields, first real, then complex, the second last the lattice site, and the first index the time step

class GPEData(object):
    def __init__(self,fName,nLat=256):
        self.latData, self.modData = readASCII_header(fName)
        self.latData['n'] = nLat
        
        self.fld = readASCII_1d(fName,nLat)
        self.nFld = self.fld.shape[-1]
        self.x = self.latData['dx']*np.arange(self.fld.shape[1])
        self.t = self.latData['dt']*np.arange(self.fld.shape[0])

        self.phi = np.angle(self.fld)
        for i in range(self.nFld):
            self.phi[:,:,i] = unwindPhase(self.phi[:,:,i])
        self.dens = np.abs(self.fld)**2

        self.dphi = np.diff(self.phi,axis=-1); self.tphi = np.mean(self.phi,axis=-1)
        return

    def getPhase(self):
        return self.phi

    def getDens(self):
        return self.dens

    def partNum(self):
        return np.sum(np.abs(self.fld)**2,axis=-1)
    
    def particleCurrent(self):
        return Particle_Current(self.fld,self.latData['dx'])

    def energy(self):
        return Energy(self.fld,self.latData['dx'],self.modData['nu'],self.modData['lam'],self.modData['om'])

    def h3(self):
        return None
    
def readASCII_header(fName):
    """
    Read header information from the files output by my code
    """
    with open(fName) as f:
        pass
    latDict = {
        "dx"    : dx,
        "dt"    : dt,
        "n"     : 256
        }
    modDict = {
        "nu"    : nu,
        "gs"    : 1.,
        "om"    : om,
        "lam"   : lv,
        "rho"   : 89.442719099991592
        }
    return latDict, modDict
    
def readASCII_1d(fName,nLat=256):
    """
    Read in simulation data stored as columns of the real and imaginary parts of the fields in alternating columns.
    
    Input :
      fName - Name of the file to read data from
      nLat - Number of lattice sites (to be removed in the future through automated method)

    Return:
      Complex array storing the value of each of the complex wavefunction condensates
    """
    a = np.loadtxt(fName,dtype=np.float64)
    nFld = a.shape[-1]/2 # Fix this when I include x and t
    return a.reshape((-1,nLat,nFld,2)).view(np.complex128)[:,:,:,0]

def readASCII_1d_real(fName,nLat=256):
    """
    Read in simulation data stored as columns of the real and imaginary parts of the fields in alternating columns, and store as a real array
    
    Input :
      fName - Name of the file to read data from
      nLat - Number of lattice sites (to be removed in the future through automated method)

    Return:
      Complex array storing the value of each of the complex wavefunction condensates
    """
    a = np.loadtxt(fName,dtype=np.float64)
    nFld = a.shape[-1]/2 # Fix this when I include x and t
    return a.reshape((-1,nLat,2*nFld))
    
def computeRho(d):
    return np.sum(np.abs(d)**2,axis=-1)

def computeRhoDiff(d):
    return np.abs(d[:,:,1])**2-np.abs(d[:,:,0])**2

def computePhi(d):
    return

def geometricDensity(d):
    """
    Returns the product of the density
    """
    return

# Check the derivative normalization in this one
def Particle_Current(d,dx):
    """
    Compute and return the particle current at each point on the lattice using complex decomposition
    """
    n = d.shape[1]; norm = 2.*np.pi/(dx)
    fk = np.fft.fft(d,axis=1)
    fk = fk*1j*np.fft.fftfreq(n)[np.newaxis,:,np.newaxis]*norm
    dd = np.fft.ifft(fk,axis=1)
    mom = np.real(d)*np.imag(dd)-np.imag(d)*np.real(dd)
    return mom

def Energy(d,dx,nu,da,om=1.,dt=1.):
    """
    Compute Energy density using complex decomposition
    """
    n = d.shape[1]; norm = 2.*np.pi/(dx)
    fk = np.fft.fft(d,axis=1)
    fk = fk*1j*np.fft.fftfreq(n)[np.newaxis,:,np.newaxis]*norm
    dd = np.fft.ifft(fk,axis=1)
    tv = np.arange(d.shape[0])*dt
    
    en_k = 0.5*np.sum(np.abs(dd)**2,axis=-1)
    en_p = 0.5*np.sum(np.abs(d)**4,axis=-1)
    en_it = -2.*om*da*(0.5*nu)**0.5*np.cos(om*tv[:,np.newaxis]))*np.real(d[:,:,0]*np.conj(d[:,:,1]))
    en_i0 = -2.*nu*np.real(d[:,:,0]*np.conj(d[:,:,1]))
    #en = 0.5*np.sum(np.abs(dd)**2,axis=-1) + 0.5*np.sum(np.abs(d)**4,axis=-1) - 2.*nu*np.real(d[:,:,0]*np.conj(d[:,:,1]))
    en = en_k + en_p + en_it + en_i0
    # nueff = nu + om*lambda*sqrt(0.5*nu)*np.cos(t)
    return en, [en_k, en_p, en_0, en_it]

def grad2_fd(d,dx):
    df2 = np.sum(np.abs(np.diff(d,axis=1))**2,axis=-1)
    dfe = np.sum(np.abs(d[:,0]-d[:,-1])**2,axis=-1)
    g2 = np.empty((d.shape[0],d.shape[1]))
    g2[:,0] = dfe + df2[:,0]
    g2[:,-1] = dfe + df2[:,-1]
    g2[:,1:-1] = df2[:,1:] + df2[:,:-1]
    return (0.5/dx**2)*g2

# Have to put in correct periodic boundary conditions
def Energy_discrete(d,dx,nu,da,om=1.,dt=1.):
    tv = dt*np.arange(d.shape[0])

    en_k = 0.5*grad2_fd(d,dx)
    en_p = 0.5*np.sum(np.abs(d)**4,axis=-1)
    en_i = -2.*(nu+om*da*(0.5*nu)**0.5*np.cos(om*tv[:,np.newaxis]))*np.real(d[:,:,0]*np.conj(d[:,:,1]))
    return en_k+en_p+en_i

def Energy_zero(d,dx,nu):
    """
    Compute zeroth order energy, droping time-dependent contribution
    """
    en = 0.5*np.sum(np.abs(dd)**2,axis=-1) + 0.5*np.sum(np.abs(d)**4,axis=-1) - 2.*nu*np.real(d[:,:,0]*np.conj(d[:,:,1]))
    return en
    
def Energy_tave(d,dx,nu,da):
    """
    Compute Effective Time-Averaged Energy Density
    """
    n = d.shape[1]; norm = 2.*np.pi/(dx)
    fk = np.fft.fft(d,axis=1)
    fk = fk*1j*np.fft.fftfreq(n)[np.newaxis,:,np.newaxis]*norm
    dd = np.fft.ifft(fk,axis=1)

    en_k = 0.5*np.sum(np.abs(dd)**2,axis=-1)
    en_p = 0.5*np.sum(np.abs(d)**4,axis=-1)
    en_i = 0.
    en = en_k + en_p + en_i
    return en

def Particle_Current_real(d,dx):
    """
    Compute and return particle current using real decomposition
    """
    n = d.shape[1]; norm = 2.*np.pi/(dx)
    fk = np.fft.rfft(d,axis=1) 
    fk = fk*1j*np.fft.rfftfreq(n)[np.newaxis,:,np.newaxis]*norm
    dd = np.fft.irfft(fk,axis=1)
    mom = [ (d[:,:,2*i]*dd[:,:,2*i+1] - d[:,:,2*i+1]*dd[:,:,2*i]) for i in range(2) ]
    return np.array(mom)

def Energy_real(d,dx):
    return

# This routine should be vectorized
from winding_number import windingNumber
def unwindPhase(f):
    """
    Unwind the phase to include the effects of winding numbers
    """
    fw = np.empty(f.shape); fw = np.mod(f+2.*np.pi,2.*np.pi)
    w0 = windingNumber(fw[:,0])
    fw[:,0] = fw[:,0]+2.*np.pi*w0  # Make left edge continuous
    w = np.empty(fw.shape)
    for i in range(len(fw)):
        w[i] = windingNumber(fw[i])
    return fw + 2.*np.pi*w

if __name__=="__main__":
    pass
#    dat = GPEData('fields.dat',1024)
#    dphi = dat.phi[:,:,1]-dat.phi[:,:,0]

    # Temporary figures
#    f,a = plt.subplots()
#    a.plot(dat.t, np.abs(np.mean(dat.energy()[0],axis=-1) - np.mean(dat.energy()[0][0])) )
#    a.set_ylabel(r'$\left|E-E_{\rm init}\right|$'); a.set_xlabel(r'$\bar{t}$')
#    f.show()

#    f,a = plt.subplots()
#    a.plot(dat.t, np.abs(np.mean(dat.partNum(),axis=-1)-np.mean(dat.partNum()[0])) )
#    a.set_ylabel(r'$\left|\rho-\rho_{\rm init}\right|$'); a.set_xlabel(r'$\bar{t}$')
#    f.show()

#    f,a = plt.subplots()
#    a.plot(dat.t, np.abs(np.mean(np.sum(dat.particleCurrent(),axis=-1),axis=-1)-np.mean(np.sum(dat.particleCurrent()[0],axis=-1),axis=-1)) )
#    a.set_ylabel(r'$\left|P-P_{\rm init}\right|$'); a.set_xlabel(r'$\bar{t}$')
#    f.show()
