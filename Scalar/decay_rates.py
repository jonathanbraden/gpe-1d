#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def decay_thresh(fName,nt,thresh=-0.8,cut=-0.5):
    a = np.genfromtxt(fName)
    d = a[:a.shape[0]-a.shape[0] % nt,:].reshape((-1,nt,5))
#    d = np.genfromtxt(fName).reshape((-1,nt,5))
    td = np.where(np.max(d[:,:,-1],axis=-1) > cut)
    t = np.argmax(d[td,:,-1] > thresh,axis=-1)[0]
    return t, d.shape[0]

def decay_thresh_interp(fName,nt,thresh=-0.8,cut=-0.5):
    """
    Determine decay times, using linear interpolation to determine first passage.  

    Returns the times in units of the output time step
    To Do: Deal with threshold index being zero better
    """
    dt = 1.

    a = np.genfromtxt(fName)
    d = a[:a.shape[0]-a.shape[0] % nt,:].reshape((-1,nt,5))
    td = np.where(np.max(d[:,:,-1],axis=-1) > cut)
    ti = np.argmax(d[td,:,-1] > thresh,axis=-1)[0]

    t = ti + dt*(cut - d[td,ti,-1]) / (d[td,ti,-1] - d[td,ti-1,-1]) # check ti isn't zero
    return t[0], d.shape[0]

def decay_thresh_notrim(fName,nt):
    d = np.genfromtxt(fName).reshape((-1,nt,5))
    t = np.argmax(d[:,:,-1] > -0.8, axis=-1)
    return t, d.shape[0]

def quick_plot(files,lv,t=-0.8,tc=-0.5):
    dt = 8.94
    t = [ decay_thresh_interp(f,100,thresh=t,cut=tc) for f in files ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax

def quick_plot_survive():
    dt = 8.94
    lv = [ '1.2','1.3','1.4','1.5','1.55' ]
#    files = [ 'bubble-count-l%s-len50.dat' % s for s in lv ]
#    files = [ 'bubble-count-l%s-len25-kc512o8.dat' % s for s in lv ]
#    files = [ 'runs/l%s/len50/bubble-count.dat' %s for s in lv ]
    files = [ 'runs/l1.2/len%s/bubble-count.dat' %s for s in ['12.5','25','50','100'] ]
    t = [ decay_thresh(f,100) for f in files ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),1.-np.linspace(0.,1.,tc[1])[:tc[0].size],label=r'$%s$' % lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax
    
def quick_plot_notrim():
    dt = 8.94
    lv = [ '1.2','1.3','1.4','1.5','1.55' ]
    files = [ 'bubble-count-l%s-len50.dat' % s for s in lv ]
    t = [ decay_thresh_notrim(f,100) for f in files ]
    for i,tc in enumerate(t):
        plt.plot(dt*np.sort(tc[0]),np.linspace(0.,1.,tc[1])[:tc[0].size],label=r'$%s$' % lv[i])
    plt.legend(loc='lower right',ncol=2)
    plt.ylabel(r'$P_{\rm decay}$')
    plt.xlabel(r'$\bar{t}$')
    
def quick_plot_len_scale(files,lv,t0=40.):
    dt = 8.94
    t = [ decay_thresh(f,100) for f in files ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot( 2**(i-3)*(dt*np.sort(tc[0])-t0),np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax

if __name__=="__main__":
    files = [ 'runs/convergence/l1.2/len12.5/kcut0.25/n2048/bubble-count.dat',
              'runs/convergence/l1.2/len12.5/kcut0.5/n2048/bubble-count.dat',
              'runs/convergence/l1.2/len12.5/kcut1/n1024/bubble-count.dat',
              'runs/convergence/l1.2/len12.5/kcut2/n1024/bubble-count.dat',
              'runs/convergence/l1.2/len12.5/kcut4/n1024/bubble-count.dat'
              ]

    base = 'runs/rough/len12.5/'
    files_r = [ base+'l%s/n1024/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.55'] ]
