#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# These functions assume that the variables are being stored in a single array with the last column indexing the fields, first real, then complex, the second last the lattice site, and the first index the time step

class GPEData(object):
    def __init__(self,fName,nLat=256):
        self.fld = readASCII_1d(fName,nLat)
        self.nFld = self.field.shape[-1]
        self.phi = np.angle(self.fld)
        self.dens = np.abs(self.fld)**2

        self.dphi = np.diff(self.phi,axis=-1); self.tphi = np.sum(self.phi,axis=-1)
        return

    def getPhase(self):
        return self.phi

    def getDens(self):
        return self.dens

    # Plotting stuff, will be nice to just define a child class
    def plotPhi(self):
        fig,ax = plt.subplots()
        return fig,ax

    def plotCosPhi(self):
        fig,ax = plt.subplots()
        return fig,ax

    def plotDens(self):
        fig,ax = plt.subplots()
        return fig,ax
    
class GPEPlot(GPEData):
    def __init__(self,fName,nLat=256):
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
    return a.reshape((-1,nLat,nFld,2)).view(np.complex128)[:,:,:,0]

def geometricDensity(a):
    """
    Returns the product of the density
    """
    return

def Particle_Current():
    """
    Compute and return the particle current at each point on the lattice
    """
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

def main():
    return

if __name__=="__main__":
    a = readASCII_1d('fields.dat',256)
    phi = np.angle(a)
    dphi = phi[:,:,1]-phi[:,:,0]
    dphi_ = unwindPhase(dphi)
