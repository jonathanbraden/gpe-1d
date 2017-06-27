#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from myplotutils import insert_rasterized_contour_plot

def readFile(fName,n):
    return np.genfromtxt(fName).reshape((-1,n,4))

def plotCosPhase(d):
    cnt = plt.contourf( (d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,1]**2+d[:,:,0]**2)*(d[:,:,2]**2+d[:,:,3]**2)),51)
    cb = plt.colorbar(label=r'$\cos(\phi)$',fraction=0.05,pad=0.01)
    plt.xlabel('position'); plt.ylabel('time')

#    for c in cnt.collections:
#        c.set_edgecolor("face")
    cb.solids.set_rasterized(True)
    return cnt
    
def plotPhase(d):
    plt.contourf( np.arccos(d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,1]**2+d[:,:,0]**2)*(d[:,:,2]**2+d[:,:,3]**2)),51)
    plt.colorbar(label=r'$\phi$',fraction=0.05,pad=0.01)
    plt.xlabel('position'); plt.ylabel('time')
    
def plotRho(d):
    plt.contourf(np.sum(d**2,axis=-1),51)
    plt.colorbar()

def plotRho1(d):
    plt.contourf(np.sum(d[:,:,2]**2,axis=-1)/np.sum(d**2,axis=-1),51)
    plt.colorbar()
    
def plotRhoTot(d):
    plt.plot(np.sum(np.sum(d**2,axis=-1),axis=-1))
    plt.xlabel('time'); plt.ylabel(r'$\rho$')
        
def main():
    return

if __name__=="__main__":
    n = 1024
    data = readFile('fields.dat',n)
    c = plotCosPhase(data)
    for _c in c.collections:
        _c.set_rasterized(True)
#        _c.set_edgecolor('Face')
