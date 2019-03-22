import gross_pitaevskii as gpe
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import palettable as pal

# Formatting parameters
cosCMap = pal.colorbrewer.diverging.RdYlBu_11_r.mpl_colormap #'coolwarm'
rhoCMap = 'PuOr'

# Labels
xLab = r'$\bar{x}$'
tLab = r'$\bar{t}$'

def plotCosPhi(d,dx,dt):
    tv = dt*np.arange(d.shape[0])
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
def noether_plot(dat,dx,dt,nu=2.e-3,da=0.,om=1.,wEn=True):
    tv = dt*np.arange(dat.shape[0])
    rho = gpe.computeRho(dat)
    mom = np.sum(gpe.Particle_Current(dat,dx),axis=-1)
    en, en_i = gpe.Energy(dat,dx,nu,da,om)

    f,a = plt.subplots()
    a.plot(tv,np.abs(np.mean(rho,axis=-1)-np.mean(rho[0])),label=r'$\varrho$')
    a.plot(tv,np.abs(np.mean(mom,axis=-1)-np.mean(mom[0])),label=r'$\mathcal{P}$')
    if wEn:
        a.plot(tv,np.mean(en,axis=-1)-np.mean(en[0]),label=r'$\mathcal{H}$')

    a.set_xlabel(tLab); a.set_ylabel(r'$\left|\mathcal{C}-\mathcal{C}_{\rm init}\right|$')
    a.set_xlim((tv[0],tv[-1]))

    f.legend(loc='upper right',bbox_to_anchor=(0,0,1,1),transform=f.transFigure)
    return f,a

def noether_plot_panel(a,dat,dx,dt,nu=2.e-3,da=0.,om=1.,wEn=True):
    tv = dt*np.arange(dat.shape[0])
    rho = gpe.computeRho(dat)
    mom = np.sum(gpe.Particle_Current(dat,dx),axis=-1)
    en, en_i = gpe.Energy(dat,dx,nu,da,om)

    a.plot(tv,np.abs(np.mean(rho,axis=-1)-np.mean(rho[0])),label=r'$\varrho$')
    a.plot(tv,np.abs(np.mean(mom,axis=-1)-np.mean(mom[0])),label=r'$\mathcal{P}$')
    if wEn:
        a.plot(tv,np.abs(np.mean(en,axis=-1)-np.mean(en[0])),label=r'$\mathcal{H}$')

    a.set_ylabel(r'$\left|\mathcal{C}-\mathcal{C}_{\rm init}\right|$')
    a.set_xlabel(tLab)
    a.set_xlim((tv[0],tv[-1]))
    return

# Write this into GPE module and then call in subroutine below
#def getPhaseDiff(dat):
#    return (d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,0]**2+d[:,:,1]**2)*(d[:,:,2]**2+d[:,:,3]**2))

def cosphi_plot_panel(a,dat,dx,dt):
    tv = dt*np.arange(dat.shape[0])
    xv = dx*np.arange(dat.shape[1])

    a.contourf(xv,tv,gpe.computeCosPhi(dat),np.linspace(-1.,1.,51),zorder=-20,cmap=cosCMap)
    a.set_rasterization_zorder(-10)

#    a.text(r'$\frac{\hbar k_{\rm nyq}}{g_{\rm s}\bar{n}} = $',0.95,0.95,va='top',ha='right',transform=transAxes)
    
    _t_x_labels(a)
    return

def rho_tot_plot_panel(a,dat,dx,dt):
    tv = dt*np.arange(dat.shape[0])
    xv = dx*np.arange(dat.shape[1])

    a.contourf(xv,tv,gpe.computeRho(dat),np.linspace(-1.,1.,51),zorder=-20,cmap=rhoCmap)
    a.set_rasterization_zorder(-10)

    _t_x_labels(a)
    return

def rho_diff_plot_panel(a,dat,dx,dt):

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
