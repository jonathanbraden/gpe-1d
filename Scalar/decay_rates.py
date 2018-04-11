#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def decay_thresh(fName,nt):
    d = np.genfromtxt(fName).reshape((-1,nt,5))
    td = np.where(np.max(d[:,:,-1],axis=-1) > -0.5)
    t = np.argmax(d[td,:,-1] > -0.8,axis=-1)[0]
    return t, d.shape[0]

def decay_thresh_notrim(fName,nt):
    d = np.genfromtxt(fName).reshape((-1,nt,5))
    t = np.argmax(d[:,:,-1] > -0.8, axis=-1)
    return t, d.shape[0]

def quick_plot():
    dt = 8.94
    lv = [ '1.2','1.3','1.4','1.5','1.55' ]
#    files = [ 'bubble-count-l%s-len50.dat' % s for s in lv ]
    files = [ 'bubble-count-l%s-len25-kc512o8.dat' % s for s in lv ]
    t = [ decay_thresh(f,100) for f in files ]
    for i,tc in enumerate(t):
        plt.plot(dt*np.sort(tc[0]),np.linspace(0.,1.,tc[1])[:tc[0].size],label=r'$%s$' % lv[i])
    plt.legend(loc='lower right',ncol=2)
    plt.ylabel(r'$P_{\rm decay}$')
    plt.xlabel(r'$\bar{t}$')
    plt.ylim(-0.1,1.1)
    
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
    
if __name__=="__main__":
    pass