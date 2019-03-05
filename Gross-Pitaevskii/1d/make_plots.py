#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from myplotutils import insert_rasterized_contour_plot
from matplotlib.ticker import MaxNLocator

# Notation.  To move to separate file
xLab = r'$\bar{x}$'
tLab = r'$\bar{t}$'

# Formatting of colours, etc.
cosCMap = 'coolwarm'
rhoCMap = 'PuOr'

lVals = [ '0','0.9','1','1.2','1.4','1.5' ]
lamFiles = [ 'paper-data/vary-lambda/fields-n512-w250-l%s.dat' % l for l in lVals ]
nLam = 512; dtLam = 1.9845101615190048; dxLam = 50./nLam
xvLam = dxLam*np.arange(nLam)

#fName = 'paper-data/vary-lambda/fields-n512-w250-l1.3.dat'
#fName = 'paper-data/fields-n256-w250-l1.4-dt1.dat'
#fName = 'spec-stencil/fields-280-rho2.dat'
#fName = 'paper-data/floquet-spec/fields-n800-w250-l1.5.dat'
#n=800; dx = 50./n; dt = 1.9845101615190048
#fName = 'paper-data/vary-lambda/fields-n512-w250-l0.dat'
#fName = 'paper-data/floquet-spec/fields-n1271-w250-l1.5-discrete.dat'
#n=490; dx = 50./n; dt = 1.9845101615190048
#n=512; dx = 50./n; dt = 1.9845101615190048
#fName = 'paper-data/floquet-spec/fields-n2048-w250-l1.5.dat'
#n=2048; dx = 50./2048; dt = 1.9845101615190048

fName0 = 'spec-stencil/fields-280-rho2.dat'
n0=280; dx0 = 1.9964892656248123; dt0 = 1.9845101615190048

fName1 = 'spec-stencil/fields-312-rho2.dat'
n1=312; dx1 = 1.7917211358171392; dt1 = 1.9845101615190048

fName2 = 'discrete-stencil/fields-440-rho2.dat'
n2=440; dx2 = 1.2704931690339714; dt2 = 1.9845101615190048

fName3 = 'discrete-stencil/fields-490-rho2.dat'
n3=490; dx3 = 1.1408510089284642; dt3 = 1.9845101615190048

n=n0; dx=dx0
xVals = dx*np.arange(n)

def multipanel_vary_lambda(dx,dt,n=512):
    """
    Make multipanel plot where we show the effect of varying lambda

    This is basically a plot as it will appear in a paper.
    """
    lVals = [ '0','0.9','1','1.2','1.4','1.5' ]
    #files = [ 'paper-data/vary-lambda/fields-n512-w250-l%s.dat' % l for l in lVals ]
    fig,ax = plt.subplots(nrows=3,ncols=2,figsize=(3.375*2.,2.08*3.),sharex=True,sharey=True,gridspec_kw = dict(wspace=0.02,hspace=0.02))

    for a,l in zip(ax.flatten(),lVals):
        f = 'paper-data/vary-lambda/fields-n512-w250-l%s.dat' % l
        d = readFile(f,n)
        tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1]) 
        cnt = a.contourf(xVals,tVals,(d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,1]**2+d[:,:,0]**2)*(d[:,:,2]**2+d[:,:,3]**2)),np.linspace(-1.,1.,51),zorder=-20,cmap=cosCMap)
        a.set_rasterization_zorder(-10)
        a.text(0.95,0.95,r'$\lambda = %s$' % l,horizontalalignment='right',verticalalignment='top',transform=a.transAxes,bbox=dict(facecolor='white'))

    # Add axis labels
    for a in ax:
        a[0].set_ylabel(tLab)
    for a in ax.T:
        a[-1].set_xlabel(xLab)

    # Fix tick marks
    for a in ax.flatten():
        a.get_xaxis().set_major_locator(MaxNLocator(5))
        a.get_yaxis().set_major_locator(MaxNLocator(5))

    # Now add the colorbar
    cb = fig.colorbar(cnt,label=r'$\cos\phi$',ticks=[-1,0,1],pad=0.05,frac=0.1)
    cb.solids.set_rasterized(True)
    return fig,ax

def multipanel_lambda_plot(dx=1.,dt=1.,n=512):
    from mpl_toolkits.axes_grid1 import ImageGrid
    fig = plt.figure(figsize=(3.375*2.,2.08*3.) )  # fix the boundary sizes
    grid = ImageGrid(fig,111,
                     nrows_ncols=(3,2),
                     axes_pad=0.05,
                     cbar_mode='single',
                     cbar_location='right',
                     cbar_pad=0.1,
                     cbar_size=0.1,
                     label_mode='L'
                     )

    for a,l in zip(grid,lVals):
        f = 'paper-data/vary-lambda/fields-n512-w250-l%s.dat' % l
        d = readFile(f,n)
        tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1]) 
        cnt = a.contourf(xVals,tVals,(d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,1]**2+d[:,:,0]**2)*(d[:,:,2]**2+d[:,:,3]**2)),np.linspace(-1.,1.,51),zorder=-20,cmap=cosCMap)
        a.set_rasterization_zorder(-10)
        a.text(0.95,0.95,r'$\lambda = %s$' % l,horizontalalignment='right',verticalalignment='top',transform=a.transAxes,bbox=dict(facecolor='white'))

    return fig, grid

def getPhase(d):
    return (d[:,:,0]*d[:,:,2] - d[:,:,1]*d[:,:,3]) / np.sqrt( (d[:,:,0]**2+d[:,:,1]**2)*(d[:,:,2]**2+d[:,:,3]**2) )

def getPhaseDiff(d):
    return (d[:,:,0]*d[:,:,2] + d[:,:,1]*d[:,:,3]) / np.sqrt( (d[:,:,0]**2+d[:,:,1]**2)*(d[:,:,2]**2+d[:,:,3]**2) )

def getRhoDiff(d):
    return d[:,:,3]**2+d[:,:,2]**2 - d[:,:,1]**2-d[:,:,0]**2

def getRhoTot(d):
    return np.sum(d**2,axis=-1)

def readFile(fName,n):
    return np.genfromtxt(fName).reshape((-1,n,4))

# Change the colour map in here
def plotCosPhase(d,dx,dt):
    f,a = plt.subplots()
    tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1])
    cnt = a.contourf(xVals, tVals, (d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,1]**2+d[:,:,0]**2)*(d[:,:,2]**2+d[:,:,3]**2)),np.linspace(-1.,1.,51),zorder=-20, cmap=cosCMap)
    #for c in cnt.collections:
    #    c.set_edgecolor("face")
    a.set_rasterization_zorder(-10)

    cb = f.colorbar(cnt,label=r'$\cos\phi$',fraction=0.05,pad=0.01,ticks=[-1,0,1])
    cb.set_label(r'$\cos\phi$',labelpad=-1)
    cb.solids.set_rasterized(True)

    _t_x_labels(a)
    return f,a

def plotPhase(d):
    f,a = plt.subplots()
    tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1])
    cnt = a.contourf(xVals, tVals, np.arccos( (d[:,:,0]*d[:,:,2]+d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,1]**2+d[:,:,0]**2)*(d[:,:,2]**2+d[:,:,3]**2)) ),51,zorder=-20)
    for c in cnt.collections:
        c.set_edgecolor("face")
    a.set_rasterization_zorder(-10)
    cb = plt.colorbar(label=r'$\phi$',fraction=0.05,pad=0.01,ticks=[0.,0.5*np.pi,np.pi])
    cb.solids.set_rasterized(True)

    _t_x_labels(a)
    return cnt

def plotTotPhase(d):
    f,a = plt.subplots()
    tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1])
    
    ph = np.arccos( (d[:,:,0]*d[:,:,2]-d[:,:,1]*d[:,:,3])/np.sqrt((d[:,:,0]**2+d[:,:,1]**2)*(d[:,:,2]**2+d[:,:,3]**2)) )
    ph_m = np.mean(ph,axis=-1)

    cnt = a.contourf(xVals,tVals,ph,51,zorder=-20)
    for c in cnt.collections:
        c.set_edgecolor("face")
    a.set_rasterization_zorder(-10)
    
    cb = plt.colorbar(cnt,label=r'$\theta - \bar{\theta}$',fraction=0.05,pad=0.01,ticks=[0.,0.5*np.pi,np.pi])
    cb.solids.set_rasterized(True)
    
    _t_x_labels(a)
    return f,a
                                                    
def plotRhoTot(d,dx,dt):
    f,a = plt.subplots()
    tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1])
    meanRho = np.mean(np.sum(d**2,axis=-1),axis=-1) 

    im = a.imshow(np.sum(d**2,axis=-1)-meanRho[:,np.newaxis],cmap=rhoCMap,extent=[xVals[0],xVals[-1],tVals[0],tVals[-1]],aspect='auto',origin='lower',vmin=-1,vmax=1)

    cb = f.colorbar(im,label=r'$\rho$',fraction=0.05,pad=0.01,ticks=[-2,-1,0,1,2])
    cb.set_label(r'$\rho-\bar{\rho}$',labelpad=-2)
    cb.solids.set_rasterized(True) 

    _t_x_labels(a)
    return f,a

def plotRhoDiff(d,dx,dt):
    f,a = plt.subplots()
    tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1])
    im = a.imshow(d[:,:,2]**2+d[:,:,3]**2-d[:,:,0]**2-d[:,:,1]**2,cmap=rhoCMap,extent=[xVals[0],xVals[-1],tVals[0],tVals[-1]],aspect='auto',origin='lower',vmin=-1,vmax=1)

    cb = f.colorbar(im,label=r'$\epsilon$',fraction=0.05,pad=0.01,ticks=[-2,-1,0,1,2])
    cb.set_label(r'$\rho_2-\rho_1$',labelpad=-2)
    cb.solids.set_rasterized(True)

    _t_x_labels(a)
    return f,a

def plotRho1(d):
    plt.contourf(np.sum(d[:,:,2]**2+d[:,:,3]**2,axis=-1)/np.sum(d**2,axis=-1),51)
    plt.colorbar()
    
def plotRhoTot_Mean(d,f=None,a=None):
    if a == None:
        f,a = plt.subplots()
    a.plot(np.sum(np.sum(d**2,axis=-1),axis=-1))
    a.set_xlabel(tLab); a.set_ylabel(r'$\rho$')
    return f,a

def plotRhoFluc_mean(d):
    f,a = plt.subplots()
    a.plot(np.std(np.sum(d**2,axis=-1),axis=-1))
    a.set_xlabel(tLab); a.set_ylabel(r'$\sigma_\rho$')

def _t_x_labels(a):
    a.set_ylabel(tLab); a.set_xlabel(xLab)
    return

def add_prelim_label(f,a):
    f.text(0.5,0.5,"PRELIMINARY",color='gray',alpha=0.9,va='center',ha='center',transform=a.transAxes,fontsize=25,fontname='Bitstream Vera Sans',weight='bold',rotation=30,style='normal',usetex=False)

def RMS_Plot(dt):
    data0 = readFile(fName0,n0)
    data1 = readFile(fName1,n1)
    f,a = plt.subplots()

    tv0 = dt0*np.arange(data0.shape[0])
    tv1 = dt1*np.arange(data1.shape[0])
    
    a.plot(tv0,np.std(np.sum(data0**2,axis=-1),axis=-1),label=r'No Floquet')
    a.plot(tv1,np.std(np.sum(data1**2,axis=-1),axis=-1),label=r'Floquet')

    a.set_xlabel(tLab); a.set_ylabel(r'$\sqrt{\delta\rho^2}$')
    a.set_xlim(tv0[0],tv0[-1])
    a.legend(loc='center right',bbox_to_anchor=(0,0,1,0.8))
    return f,a

def mean_cos_phi_plot(fn1,n1,fn2,n2,fn3,n3,dt):
    data1 = readFile(fn1,n1)
    data2 = readFile(fn2,n2)
    data3 = readFile(fn3,n3)

    tv = dt*np.arange(data1.shape[1])
    f1 = getPhaseDiff(data1)
    f2 = getPhaseDiff(data2)
    f3 = getPhaseDiff(data3)

    f,a = plt.subplots()
    a.plot(tv,np.mean(f1,axis=-1),label=r'$dx = $')
    a.plot(tv,np.mean(f2,axis=-1),label=r'$dx = $')
    a.plot(tv,np.mean(f3,axis=-1),label=r'$dx = $')

    return f,a
    
def relative_stats_plot():
    data0 = readFile(fName0,n0)
    f0 = getPhaseDiff(data0); tv0 = dt0*np.arange(data0.shape[0])
    data1 = readFile(fName1,n1)
    f1 = getPhaseDiff(data1); tv1 = dt1*np.arange(data1.shape[0])

    # Mean cos of relative phase
    f,a = plt.subplots()
    a.plot(tv0,np.mean(f0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.mean(f1,axis=-1),label=r'Floquet')
    a.set_xlabel(tLab); a.set_ylabel(r'$\langle\cos\phi\rangle_{\rm V}$')

    a.set_xlim((tv0[0],tv0[-1])); a.set_ylim((-1.1,1.1))
    a.legend(loc='lower right')
    f.savefig('mean-cos-phi.pdf')
    f.show()
    
    # std dev of cos of relative phase
    f,a = plt.subplots()
#    a.plot(tv0,np.mean(f0,axis=-1),label=r'No Floquet')
#    a.plot(tv1,np.mean(f1,axis=-1),label=r'Floquet')
    a.plot(tv0,np.std(f0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.std(f1,axis=-1),label=r'Floquet')
    a.set_xlabel(tLab); a.set_ylabel(r'$\langle(\delta\cos\phi)^2\rangle_{\rm V}$')

    a.set_xlim((tv0[0],tv0[-1])); a.set_ylim((-1.1,1.1))
    a.legend(loc='lower right')
    f.savefig('std-cos-phi.pdf')
    f.show()

    f,a = plt.subplots()
    e0 = getRhoDiff(data0); e1 = getRhoDiff(data1)
    a.plot(tv0,np.mean(e0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.mean(e1,axis=-1),label=r'Floquet')
    a.set_xlabel(tLab); a.set_ylabel(r'$\langle\epsilon\rangle_{\rm V}$')
    a.set_xlim((tv0[0],tv0[-1]))
    f.savefig('mean-rho-diff.pdf')
    f.show()

    f,a = plt.subplots()
    a.plot(tv0,np.std(e0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.std(e1,axis=-1),label=r'Floquet')
    a.set_xlabel(tLab); a.set_ylabel(r'$\langle\delta\epsilon^2\rangle_{\rm V}$')
    a.set_xlim((tv0[0],tv0[-1]))
    a.legend(loc='center right')
    f.savefig('std-rho-diff.pdf')
    f.show()

    f,a = plt.subplots()
    r0 = getRhoTot(data0); r1 = getRhoTot(data1)
    a.plot(tv0,np.std(r0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.std(r1,axis=-1),label=r'Floquet')
    a.set_xlabel(tLab); a.set_ylabel(r'$\langle\delta\rho^2\rangle_{\rm V}$')
    a.set_xlim((tv0[0],tv0[-1]))
    a.legend(loc='center right')
    f.savefig('std-rho-tot.pdf')
    f.show()

    f,a = plt.subplots()
    a.plot(tv0,np.std(e0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.std(e1,axis=-1),label=r'Floquet')
    a.plot(tv0,np.std(r0,axis=-1),label=r'No Floquet')
    a.plot(tv1,np.std(r1,axis=-1),label=r'Floquet')
    a.set_xlabel(tLab); a.set_ylabel(r'$\langle\delta\rho^2\rangle_{\rm V}$')
    a.set_xlim((tv0[0],tv0[-1]))
    a.legend(loc='center right')
    f.savefig('std-rho-both.pdf')
    f.show()
    
    return
    
def floquet_plots():
    data0 = readFile(fName0,n0)
    f,a = plotCosPhase(data0,dx0,dt0)
    f.savefig('bubbles-no-floquet-spec.pdf',dpi=300)

    data1 = readFile(fName1,n1)
    f,a = plotCosPhase(data1,dx1,dt1)
    f.savefig('bubbles-w-floquet-spec.pdf',dpi=300)
    
    data2 = readFile(fName2,n2)
    f,a = plotCosPhase(data2,dx2,dt2)
    f.savefig('bubbles-no-floquet-discrete.pdf',dpi=300)
    
    data3 = readFile(fName3,n3)
    f,a = plotCosPhase(data3,dx3,dt3)
    f.savefig('bubbles-w-floquet-discrete.pdf',dpi=300)
    
    return

def getParticleSpectrum(fName,n):
    # Fix this so I'm not getting the reflected spectrum
    data = readFile(fName,n)
    
    f1 = data[:,:,0]+1j*data[:,:,1]; f2 = data[:,:,2]+1j*data[:,:,3]
    fk1 = np.fft.fft(f1,axis=-1)
    fk2 = np.fft.fft(f2,axis=-1)
    
    return fk1, fk2

def plotFourierGrowth():
    data = readFile(fName1,n1)
    drho = data[:,:,2]**2+data[:,:,3]**2-data[:,:,1]**2-data[:,:,0]**2
    rho = np.sum(data**2,axis=-1)
    tv = dt1*np.arange(data.shape[0])
    
    fk_d = np.fft.rfft(drho)
    fk_t = np.fft.rfft(rho)

    f,a = plt.subplots()
    kind = [50,149,150]
    for kc in kind:
        a.plot(tv,np.abs(fk_d[:,kc]**2))
    a.set_xlabel(tLab); a.set_ylabel(r'$\left|\delta\tilde{\epsilon}_{\rm k}\right|^2$')
    a.set_yscale('log')
    a.set_xlim((0.,150.)); a.set_ylim((1.e-1,1.e3))
    f.show()
    f.savefig('fourier-modes-floquet-rel.pdf')

    f,a = plt.subplots()
    kind = [50,149,150]
    for kc in kind:
        a.plot(tv,np.abs(fk_t[:,kc]**2))
    a.set_xlabel(tLab); a.set_ylabel(r'$\left|\delta\tilde{\rho}_{\rm k}\right|^2$')
    a.set_yscale('log');
    a.set_xlim((0.,150.)); a.set_ylim((1.e-1,1.e3))
    f.show()
    f.savefig('fourier-modes-floquet-tot.pdf')

    return

def plot_variance_w_density():
    files = [ 'fields-rho1.dat', 'fields-rho100.dat', 'fields-rho10000.dat' ]
    data = np.array([ readFile(f,1024) for f in files])
    return

def main():
    return

if __name__=="__main__":
#    fNameL = 'paper-data/vary-lambda/fields-n512-w250-l0.dat'
#    nL = 512
#    dtL = 1.9845101615190048
#    dxL = 1.0918300671385692

    n = 290
    data = readFile('fields.dat',1024)
#    data = readFile('paper-data-new/convergence/floquet-dt/fields-floquet-n290-64wosc.dat',n)
    f,a = plotCosPhase(data,1.9276448081894739,1.4049629462081452)
    rho = getRhoTot(data)
    drho = getRhoDiff(data)
    phi = getPhaseDiff(data)
    
#    data = readFile('fields-n512-w256.dat',512)
#    f,a = plotCosPhase(data,1.091830067,0.175620368276)

#    f,a = RMS_Plot(dt)
    
#    f,a = plotCosPhase(data)
#    f,a = plotRhoTot(data)
#    c = plotCosPhase(data)
    pass
