#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

import gross_pitaevskii as gpe

from matplotlib.ticker import MaxNLocator

from notation import *
import palettable as pal
from plot_sizes import add_cb_above

linewidth = 5.93
halfwidth = 2.97

l_base=50./np.sqrt(2.e-3)/2.
# Compute the dx in here
dict_floq_ir_old = dict(fn='spec-stencil/fields-280-rho2.dat',
                    n=280,
                    dx=l_base/280,
                    dt=1.9845101615190048)
dict_floq_uv_old = dict(fn='spec-stencil/fields-312-rho2.dat',
                    n=312,
                    dx=l_base/312,
                    dt=1.9845101615190048)
dict_nofloq_l0 = dict(fn='paper-data/vary-lambda/fields-n512-w250-l0.dat',
                      n=512,
                      dx=l_base/512,
                      dt = 1.9845101615190048)
dict_discrete_ir = dict(fn='discrete-stencil/fields-440-rho2.dat',
                        n=440,
                        dx=l_base/440,
                        dt=1.9845101615190048)
dict_discrete_uv = dict(fn='discrete-stencil/fields-490-rho2.dat',
                        n=490,
                        dx=l_base/490,
                        dt=1.9845101615190048)

dict_floq_ir_new = dict(fn='paper-data-new/convergence/floquet-dt/fields-floquet-n295-128wosc.dat',
                        n=295,
                        dx=l_base/295,
                        dt=1.4049629462081452)

dict_floq_uv_new = dict(fn='paper-data-new/convergence/floquet-dt/fields-floquet-n312-128wosc.dat',
                        n=312,
                        dx=l_base/312,
                        dt=1.4049629642081452)


dict_floq_uv = dict_floq_uv_old
dict_floq_ir = dict_floq_ir_old

# Temporary version, copied, needs massaging
def multipanel_vary_lambda(dx=1.0918300671385692,dt=1.9845101615190048,n=512):
    """
    Make multipanel plot where we show the effect of varying lambda
    """
    lVals = [ '0','0.9','1','1.2','1.4','1.5' ]
    # files = [ 'paper-data/vary-lambda/fields-n512-w250-l%s.dat' % l for l in lVals ]

    #fig,ax = plt.subplots(nrows=3,ncols=2,sharex=True,sharey=True,gridspec_kw=dict(wspace=0.04,hspace=0.04))
    fig,ax = plt.subplots(nrows=3,ncols=2,figsize=(5.93,5.0811185271229),sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.04))
    fig.subplots_adjust(left=0.12062019861345326,right=0.9097301854974705,bottom=0.091569951634494468,top=0.94533137215851082)

    for a,l in zip(ax.flatten(),lVals):
        f = 'paper-data/vary-lambda/fields-n512-w250-l%s.dat' % l
        d = gpe.readASCII_1d(f,n)
        tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1]) 
        cnt = a.contourf(xVals,tVals,gpe.computeCosPhi(d),np.linspace(-1.,1.,51),zorder=-20,cmap=cosCMap)
        a.set_rasterization_zorder(-10)
        a.text(0.95,0.95,r'$%s = %s$' % (lPar,l),horizontalalignment='right',verticalalignment='top',transform=a.transAxes,bbox=dict(facecolor='white'))

    # Add axis labels
    for a in ax:
        a[0].set_ylabel(tLab)
    for a in ax.T:
        a[-1].set_xlabel(xLab)

    # Fix tick marks
    for a in ax.flatten():
        #a.get_xaxis().set_major_locator(MaxNLocator(5))
        #a.get_yaxis().set_major_locator(MaxNLocator(4))
        a.xaxis.set_ticks([0,100,200,300,400,500])
        a.yaxis.set_ticks([0,100,200,300])
        
    # Add colorbar
    cax = add_cb_above(fig,0.1,0.05,nr=3,nc=2,wspace=0.04,hspace=0.04)
    cb = fig.colorbar(cnt,cax=cax,orientation='horizontal',ticks=[-1,0,1],label=r'$\cos%s$' % phaseRel)
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.xaxis.set_label_position('top')
    cb.solids.set_rasterized(True)
    return fig,ax

# Should reproduce the above
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
        d = gpe.readASCII_1d(f,n)
        tVals = dt*np.arange(d.shape[0]); xVals = dx*np.arange(d.shape[1]) 
        cnt = a.contourf(xVals,tVals,gpe.computeCosPhi(d),np.linspace(-1.,1.,51),zorder=-20,cmap=cosCMap)
        a.set_rasterization_zorder(-10)
        a.text(0.95,0.95,r'$%s = %s$' % (lPar,l),horizontalalignment='right',verticalalignment='top',transform=a.transAxes,bbox=dict(facecolor='white'))
    return fig, grid

def multipanel_phase_plot(d0,d1):
    dicts = [d0,d1]
    
    f,a = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.02))

    for i in range(2):
        a_ = a[i]; d_ = dicts[i]
        dx_ = d_['dx']; dt_ = d_['dt']
        data = gpe.readASCII_1d(d_['fn'],d_['n'])
        tv = dt_*np.arange(data.shape[0]); xv = dx_*np.arange(data.shape[1])
        cnt = a_.contourf(xv,tv,gpe.computeCosPhi(data),np.linspace(-1.,1.,51),zorder=-20,cmap=cosCMap)
        a_.set_rasterization_zorder(-10)
        a_.text(0.95,0.95,r'$%s = %.2f$' % (kDim,np.pi/dx_),va='top',ha='right',transform=a_.transAxes,bbox=dict(facecolor='white',alpha=0.9))

    _2_panel_labels(a)
    # Add colorbar
    cax = add_cb_above(f,0.1,0.05)
    cb = f.colorbar(cnt,cax=cax,orientation='horizontal',ticks=[-1,0,1],label=r'$\cos%s$' % phaseRel)
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.xaxis.set_label_position('top')
    return f,a

def multipanel_phase_plot_spectral():
    f,a = multipanel_phase_plot(d0=dict_floq_ir,d1=dict_floq_uv)
    return f,a
    
def multipanel_phase_plot_discrete():
    f,a = multipanel_phase_plot(d0=dict_discrete_ir,d1=dict_discrete_uv)
    return f,a
    
def multipanel_rho_plot_no_floquet():
    d_ = dict_nofloq_l0
    dx = d_['dx']; dt = d_['dt']
    
    data = gpe.readASCII_1d(d_['fn'],d_['n'])
    tv = dt*np.arange(data.shape[0]); xv = dx*np.arange(data.shape[1])
     
    f,a = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.02))

    drho = gpe.computeRhoDiff(data)
    im = a[0].imshow(drho-np.mean(drho,axis=-1)[:,np.newaxis],cmap=rhoCMap,extent=[xv[0],xv[-1],tv[0],tv[-1]],aspect='auto',origin='lower',vmin=-1,vmax=1)
    a[0].set_rasterization_zorder(-10)
    
    rho = gpe.computeRho(data)
    im= a[1].imshow(rho-np.mean(rho,axis=-1)[:,np.newaxis],cmap=rhoCMap,extent=[xv[0],xv[-1],tv[0],tv[-1]],aspect='auto',origin='lower',vmin=-1,vmax=1)
    a[1].set_rasterization_zorder(-10)

    a[0].text(0.95,0.95,r'$%s-\langle%s\rangle_{\rm V}$' % (densRel,densRel),va='top',ha='right',transform=a[0].transAxes,bbox=dict(facecolor='white',alpha=0.9))
    a[1].text(0.95,0.95,r'$2\left(%s-\langle%s\rangle_{\rm V}\right)$' % (densTot,densTot),va='top',ha='right',transform=a[1].transAxes,bbox=dict(facecolor='white',alpha=0.9))

    _2_panel_labels(a)

    # Add colorbarx
    cax = add_cb_above(f,0.1,0.05)
    cb = f.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,0,1],label=r'$2%s,%s$' % (densTot,densRel))
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.xaxis.set_label_position('top')
    cb.solids.set_rasterized(True)
    return f,a
    
def multipanel_rho_diff_plot(d0,d1):
    dicts = [d0,d1]
    
    f,a = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.02))

    for i in range(2):
        a_=a[i]; d_=dicts[i]
        data = gpe.readASCII_1d(d_['fn'],d_['n'])
        dt_ = d_['dt']; dx_ = d_['dx']
        tv = dt_*np.arange(data.shape[0]); xv = dx_*np.arange(data.shape[1])
        drho = gpe.computeRhoDiff(data)
        im = a_.imshow(drho-np.mean(drho,axis=-1)[:,np.newaxis],cmap=rhoCMap,extent=[xv[0],xv[-1],tv[0],tv[-1]],aspect='auto',origin='lower',vmin=-1,vmax=1)
        a_.set_rasterization_zorder(-10)
        a_.text(0.95,0.95,r'$%s = %.2f$' % (kDim,np.pi/dx_),va='top',ha='right',transform=a_.transAxes,bbox=dict(facecolor='white',alpha=0.9))

    _2_panel_labels(a)
    # Add a colorbar
    cax = add_cb_above(f,0.1,0.05)
    cb = f.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,0,1],label=r'$%s-\langle%s\rangle_{\rm V}$' % (densRel,densRel))
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.xaxis.set_label_position('top')
    cb.solids.set_rasterized(True)
    return f,a

def multipanel_rho_diff_plot_spectral():
    f,a = multipanel_rho_diff_plot(dict_floq_ir,dict_floq_uv)
    return f,a
    
def multipanel_rho_tot_plot_spectral():
    f,a = multipanel_rho_tot_plot(dict_floq_ir,dict_floq_uv)
    return f,a

def multipanel_rho_tot_plot(d0,d1):
    dicts = [d0,d1]
    
    f,a = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.02))
        
    for i in range(2):
        a_=a[i];d_=dicts[i]
        data = gpe.readASCII_1d(d_['fn'],d_['n'])
        dt_=d_['dt']; dx_=d_['dx']
        tv = dt_*np.arange(data.shape[0]); xv = dx_*np.arange(data.shape[1])
        rho = gpe.computeRho(data)
        im = a_.imshow(rho-np.mean(rho,axis=-1)[:,np.newaxis],cmap=rhoCMap,extent=[xv[0],xv[-1],tv[0],tv[-1]],aspect='auto',origin='lower',vmin=-1,vmax=1)
        a_.set_rasterization_zorder(-10)
        a_.text(0.95,0.95,r'$%s = %.2f$' % (kDim,np.pi/dx_),va='top',ha='right',transform=a_.transAxes,bbox=dict(facecolor='white',alpha=0.9))

    _2_panel_labels(a)    
    # Add colorbar
    cax = add_cb_above(f,0.1,0.05)
    cb = f.colorbar(im,cax=cax,orientation='horizontal',ticks=[-1,0,1],label=r'$2\left(%s-\langle%s\rangle_{\rm V}\right)$' % (densTot,densTot))
    cb.ax.xaxis.set_ticks_position('bottom')
    cb.ax.xaxis.set_label_position('top')
    cb.solids.set_rasterized(True)
    
    return f,a

def _2_panel_labels(a):
    a[0].set_xlabel(xLab); a[0].set_ylabel(tLab)
    a[1].set_xlabel(xLab)
    return

# Temporary but nearing complete version of this plot
def RMS_plot_floquet():
    d0 = dict_floq_ir; d1 = dict_floq_uv
    dicts = [d0,d1]

    f,a = plt.subplots()
    
    data0 = gpe.readASCII_1d(d0['fn'],d0['n']); data1 = gpe.readASCII_1d(d1['fn'],d1['n'])
    tv0 = d0['dt']*np.arange(data0.shape[0]); tv1 = d1['dt']*np.arange(data1.shape[0])
   
    cmap = plt.get_cmap("tab10")
    # add a loop
    a.plot(tv0,np.std(gpe.computeRhoDiff(data0),axis=-1),'-',color=cmap(0),label=r'$%.2f, %s$' % (d0['dx'],densRel))
    a.plot(tv1,np.std(gpe.computeRhoDiff(data1),axis=-1),'--',lw=1.,color=cmap(2),label=r'$%.2f, %s$' % (d1['dx'],densRel) )
    a.plot(tv0,np.std(gpe.computeRho(data0),axis=-1),'-',color=cmap(1),label=r'$%.2f, %s$' % (d0['dx'],densTot) )
    a.plot(tv1,np.std(gpe.computeRho(data1),axis=-1),'--',color=cmap(3),lw=1.,label=r'$%.2f, %s$' % (d1['dx'],densTot) )

    a.set_xlabel(tLab);
    a.set_ylabel(r'$\sigma_{%s},2\sigma_{%s}$' % (densRel,densTot))
    a.set_xlim(tv0[0],tv0[-1]); a.set_ylim(-0.01,0.6)
    a.xaxis.set_major_locator(MaxNLocator(4))
    a.legend(loc='upper left',bbox_to_anchor=(0,0,1,1.18),title=r'$%s\, dx$' % xDim,ncol=2,columnspacing=1.,handlelength=1.5)
    return f,a

def RMS_plot_no_floquet():
    fName = 'paper-data/vary-lambda/fields-n512-w250-l0.dat'
    n = 512; dt = 1.9845101615190048; dx = 1.0918300671385692

    data = gpe.readASCII_1d(fName,n)
    tv = dt*np.arange(data.shape[0])

    f,a = plt.subplots()
    a.plot(tv,np.std(gpe.computeRhoDiff(data),axis=-1),label=r'$%s$' % densRel)
    a.plot(tv,np.std(gpe.computeRho(data),axis=-1),'--',label=r'$%s$' % densTot)
    
    a.set_xlabel(tLab)
    a.set_ylabel(r'$\sigma_{%s},2\sigma_{%s}$' % (densRel,densTot))
    a.set_xlim(tv[0],tv[-1])
    a.xaxis.set_major_locator(MaxNLocator(4))
    a.set_ylim((-0.005,0.2)); a.set_yticks([0.,0.1,0.2])
    a.legend(loc='lower center',bbox_to_anchor=(0,0,1,1),ncol=2)
    return f,a

def mean_cos_phi_plot():
    d0 = dict_floq_ir; d1 = dict_floq_uv
#    d0 = dict_floq_ir_old; d1 = dict_floq_uv_old
    dicts = [d0,d1]
    data1 = gpe.readASCII_1d(d0['fn'],d0['n'])
    data2 = gpe.readASCII_1d(d1['fn'],d1['n'])

    tv = d1['dt']*np.arange(data1.shape[0])
    f1 = gpe.computeCosPhi(data1)
    f2 = gpe.computeCosPhi(data2)

    f,a = plt.subplots()
    a.plot(tv,np.mean(f1,axis=-1),'-',label=r'$%.2f\ {\rm (No\ Floquet)}$' % d0['dx'])
    a.plot(tv,np.mean(f2,axis=-1),'--',label=r'$%.2f\ {\rm (Floquet)}$' % d1['dx'])
    
    a.legend(loc='lower right',bbox_to_anchor=(0,0,1,1),title=r'$%s\, dx$' % xDim)
    a.set_ylabel(r'$\left\langle \cos %s \right\rangle_{\rm V}$' % phaseRel)
    a.set_xlabel(tLab)

    a.set_yticks([-1,0,1]); a.set_ylim(-1.2,1.2)
    a.set_xlim((tv[0],tv[-1]))
    a.xaxis.set_major_locator(MaxNLocator(4))

    # Inline the phi_fv, phi_tv and phi_random lines
    a.axhline(1.,color='grey',linestyle='-.',zorder=-5,alpha=0.5)
    a.axhline(-1.,color='grey',linestyle='-.',zorder=-5,alpha=0.5)
    a.axhline(0.,color='grey',linestyle='-.',zorder=-5,alpha=0.5)

    bb = dict(pad=4.,facecolor='white',alpha=1.,edgecolor='white')
    a.text(tv[tv.size/4],1.,r'$%s_{\rm tv}$' % phaseRel,va='center',bbox=bb,zorder=-4)
    a.text(tv[tv.size/4],0.,r'random',va='center',bbox=bb,zorder=-4)
    a.text(tv[tv.size/4],-1.,r'$%s_{\rm fv}$' % phaseRel,va='center',bbox=bb,zorder=-4)
    
    return f,a

from convergence_testing import *
# Need to correct the colour cycle here
def temporal_convergence_plot():
    f,a = plt.subplots(ncols=2,nrows=1,sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.02))

    dat = collectDataFixedN(files_no_floquet_varydt_512,512)
    dt = 1.4049629462081452  # Fix this
    plot_t_convergence_axis(a[0],dat,dt)
    
    dat = collectDataFixedN(files_floquet_varydt_1024,1024)  # the 2048 doesn't read properly for some reason.  WTF!!!  Probably they are different sizes or something ...
    dt = 1.4049629462081452  # Fix this
    plot_t_convergence_axis(a[1],dat,dt)

    for a_ in a:
        a_.set_xlabel(tLab)
    a[0].set_ylabel(r'$\left|\psi^{\rm (n)}-\psi^{\rm (n-1)}\right|_{\rm max}$')

    a[0].set_ylim(1.e-15,10.)
    return f,a

def spatial_convergence_plot():
    f,a = plt.subplots(ncols=2,nrows=1,sharex=True,sharey=True,gridspec_kw = dict(wspace=0.04,hspace=0.02))

    dat = collectDataFixedT(files_no_floquet_varydx_295,295)
    dt = 1.4049629462081452  # Fix this
    plot_x_convergence_axis(a[0],dat,dt)
    
    dat = collectDataFixedT(files_floquet_varydx_295_w256,295)
    dt = 1.4049629462081452  # Fix this
    plot_x_convergence_axis(a[1],dat,dt)

    for a_ in a:
        a_.set_xlabel(tLab)
    a[0].set_ylabel(r'$\left|\psi^{\rm (n)}-\psi^{\rm (n-1)}\right|_{\rm max}$')

    a[0].set_ylim(1.e-15,10.)
    
    return f,a

def make_paper_figures():
    # Figures full width with 2 panels and no colorbar
    # Add automatic setting of figure size in here
    f,a = spatial_convergence_plot()
    f.savefig('gpe-convergence-dx.pdf')

    f,a = temporal_convergence_plot()
    f.savefig('gpe-convergence-dt.pdf')

    # Figures full width with 2 panels and a colorbar
    f,a = multipanel_phase_plot_discrete()
    f.savefig('bubbles-discrete-show-floquet.pdf')

    f,a = multipanel_phase_plot_spectral()
    f.savefig('bubbles-spectral-show-floquet.pdf')

    f,a = multipanel_rho_diff_plot_spectral()    
    f.savefig('rho-diff-evolution-floquet.pdf')

    f,a = multipanel_rho_tot_plot_spectral()
    f.savefig('rho-tot-evolution-floquet.pdf')

    f,a = multipanel_rho_plot_no_floquet()
    f.savefig('rho-evolution-l0.pdf')
    
    # Single figures 3/4 width
    # Add automatic sizing of figure in here
#    f.savefig('standard-rho-both.pdf')
    
#    f.savefig('density-rms-l0.pdf')

#    f.savefig('mean-cos-phi.pdf')

#    f.savefig('floquet-chart-fiducial.pdf')

#    f.savefig('effective-potentials.pdf')
    
    # The large 2x3 figure with a colorbar
    f,a = multipanel_vary_lambda(dx=1.0918300671385692,dt=1.9845101615190048)
    f.savefig('gpe-varyl-multipanel.pdf')
    
    return
    
if __name__=="__main__":
#    f,a = mean_cos_phi_plot()
#    f.savefig('mean-cos-phi.pdf')

#    f,a = RMS_plot_floquet()
    pass
