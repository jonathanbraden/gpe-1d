import gross_pitaevskii as gpe
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# Formatting parameters
cosCMap = 'coolwarm'
rhoCMap = 'PuOr'

# Labels
xLab = r'$\bar{x}$'
tLab = r'$\bar{t}$'

def plotCosPhi(d):
    return

def plotPhiDiff(d):
    return

def plotPhiTot(d):
    return

def plotRhoDiff(d):
    return

def plotRhoTot(d):
    return

# Add model parameters for energy
def noether_plot(dat,dx,dt,nu=2.e-3,da=0.,om=1.):
    tv = dt*np.arange(dat.shape[0])
    rho = gpe.computeRho(dat)
    mom = np.sum(gpe.Particle_Current(dat,dx),axis=-1)
    en, en_i = gpe.Energy(dat,dx,nu,da,om)

    f,a = plt.subplots()
    a.plot(tv,np.abs(np.mean(rho,axis=-1)-np.mean(rho[0])),label=r'$\rho$')
    a.plot(tv,np.abs(np.mean(mom,axis=-1)-np.mean(mom[0])),label=r'$P$')
    a.plot(tv,np.mean(en,axis=-1)-np.mean(en[0]),label=r'$H$')

    a.set_xlabel(r'$\bar{t}$'); a.set_ylabel(r'$\left|\mathcal{C}-\mathcal{C}_{\rm init}\right|$')
    a.set_xlim((tv[0],tv[-1]))
    return f,a
    

# Utility subroutines to minimize code duplication

def _t_x_labels(a):
    a.set_xlabel(xLab); a.set_ylabel(tLab)
    return

if __name__ == "__main__":
    f_no = 'paper-data-new/convergence/no-floquet-dx/fields-floquet-n4096-wosc512.dat'
    f_yes_4720 = 'paper-data-new/convergence/floquet-dx/fields-floquet-n4720-w256.dat'
    f_yes_2360 = 'paper-data-new/convergence/floquet-dx/fields-floquet-n2360-w256.dat'
    f_yes_1180_256 = 'paper-data-new/convergence/floquet-dx/fields-floquet-n1180-w256.dat'
    f_yes_1180_128 = 'paper-data-new/convergence/floquet-dx/fields-floquet-n1180-w128.dat'
    
    #n = 4096; dat = gpe.readASCII_1d(f_no,4096)
    n = 4720; dat = gpe.readASCII_1d(f_yes_4720,n)

    #dx=0.13647875839232115
    nu = 2.e-3
    dx = 50./n/2./np.sqrt(nu)  
    dt= 1.4049629462981452
    da = 1.3 #0.
    w = 50.*2.*np.sqrt(nu)

    f,a = noether_plot(dat,dx,dt,nu,da,w)
    pass
