#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# These functions assume that the variables are being stored in a single array with the last column indexing the fields, first real, then complex, the second last the lattice site, and the first index the time step

def main():
    return

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
    return a.reshape((-1,nLat,nFld,2)).view(np.complex128)

def Local_Density(a):
    return np.sum(a.view(np.float64)**2,axis=-1)

def Total_Density(a):
    """
    Return the total density in the BEC.  If the GPE is being satisfied this must be constant (at least in the homogeneous limit).  Is it still conserved in the inhomogeneous limit?"
    """
    return np.sum(Local_Density(a),axis=-1) / np.float64(a.shape[1])

def Particle_Current():
    """
    Compute and return the particle current at each point on the lattice
    """
    return

def Complex_to_Phase(psi):
    """
    Given the cartesian representations of the the BEC condensates, convert to the number phase representation
      \psi_i = \sqrt{\rho_i}e^{i\phi_i}
    The returned phase is adjusted to be continuous in time and $\rho_i$ is assumed positive.
    """
    return np.abs(psi), np.angle(psi)

def Phase_to_Cartestion(rho,phi):
    """
    Convert the number-phase representation of the BEC into a real-complex representation
    """
    return rho*(np.cos(phi) + 1j*np.sin(phi))

def Unwrap_Phase(phi):
    """
    Unwrap a stream of phase variables to be continous starting from the first index.
    """
    # Allow for phase going from pi to -pi
    alpha = where(np.abs(np.diff(phi) > 0.5*np.pi))
    # Allow for phase going from -pi to pi
    return

def Phase_Window(phi,phiMin):
    """
    Move phases to within the window [phiMin,phiMin+2*pi).
    Due to roundoff errors that are not checked for, this method may fail in some
    exceptional circumstances

    Currently will only worky within a window plus or minus 2pi of the bounds
    """
    twopi = 2.*np.pi; twopi_i = 1./twopi
    phiMax = phiMin + twopi
    wind = np.array((phi-phiMin)*twopi_i,dtype=np.int)
    phi_ = np.copy(phi)

    mask = np.where(phi > phiMax)
    phi_[mask] -= 2.*np.pi*wind[mask]
    mask = np.where(phi < phiMin)
    phi_[mask] += 2.*np.pi*(wind[mask]+1)
    return phi_

# I can combine the two unwrap functions above by simply passing in a time vs a space slice

def Relative_Phase():
    return

def Phase_to_Cartesian():
    return

if __name__=="__main__":
    a = readASCII_1d('fields.dat')
    phi = np.angle(a)[:,:,:,0]
    dphi = phi[:,:,1]-phi[:,:,0]
    rho = Local_Density(a)
    rho_0 = np.mean(np.sum(rho[0,...],axis=-1))

    plt.plot(np.abs(np.mean(np.sum(rho,axis=-1),axis=-1)-rho_0))
    a = readASCII_1d('fields_0.dat')
    rho = Local_Density(a)
    rho_0 = np.mean(np.sum(rho[0,...],axis=-1))
    plt.plot(np.abs(np.mean(np.sum(rho,axis=-1),axis=-1)-rho_0))

    plt.yscale('log')
