#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def decay_thresh(fName,nt,thresh=-0.8,cut=-0.5):
    a = np.genfromtxt(fName)
    d = a[:a.shape[0]-a.shape[0] % nt,:].reshape((-1,nt,5))
    td = np.where(np.max(d[:,:,-1],axis=-1) > cut)
    t = np.argmax(d[td,:,-1] > thresh,axis=-1)[0]
    return t, d.shape[0]

def decay_thresh_interp(fName,nt,thresh=-0.8,cut=-0.5,dt=1.):
    """
    Determine decay times, using linear interpolation to determine first passage.  

    Returns the times in units of the output time step
    To Do: Deal with threshold index being zero better
    """
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

def get_trajectories(files,nt,ax=-1):
    d = []
    for f in files:
        a = np.genfromtxt(f)
        d.append( a[:a.shape[0]-a.shape[0]% nt,ax].reshape((-1,nt)) )
    return d

def quick_plot_traj(d,lv,t=-0.8,tc=-0.5,dt=1.5440808887540916*2.*(2.e-3)**0.5,interp=True):
    t = [ extract_decay_times(dc,th=t,cut=tc,interp=interp) for dc in d ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax

def quick_plot_traj_survive(d,lv,t=-0.8,tc=-0.5,dt=1.5440808887540916*2.*(2.e-3)**0.5,interp=True):
    t = [ extract_decay_times(dc,th=t,cut=tc,interp=interp) for dc in d ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),1.-np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm survive}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax

def quick_plot(files,lv,nt=256,t=-0.8,tc=-0.5,dt = 1.5440808887540916*2.*(2.e-3)**0.5,interp=True):
    if interp:
        t = [ decay_thresh_interp(f,nt,thresh=t,cut=tc) for f in files ]
    else:
        t = [ decay_thresh(f,nt,thresh=t,cut=tc) for f in files ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    _decay_labels(ax)
    ax.set_ylim(-0.1,1.1)

    return fig,ax

def quick_plot_survive(files,lv,nt=512,t=-0.8,tc=-0.5,dt=1.5440808887540916/2./(2.e-3)**0.5,interp=True):
    if interp:
        t = [ decay_thresh_interp(f,nt,thresh=t,cut=tc) for f in files ]
    else:
        t = [ decay_thresh(f,nt,thresh=t,cut=tc) for f in files ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),1.-np.linspace(0.,1.,tc[1])[:tc[0].size],label=r'$%s$' % lv[i])
    plt.legend(loc='lower right',ncol=2)
    _survive_labels(ax)
    ax.set_ylim(-0.1,1.1)
    return fig,ax

def _decay_labels(ax):
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    return

def _survive_labels(ax):
    ax.set_ylabel(r'$P_{\rm survive}$')
    ax.set_xlabel(r'$\bar{t}$')
    return
                 
def decay_ind(d,cut):
    return np.where(np.max(d[:,:],axis=-1) > cut)

def extract_decay_times(d,th=-0.8,cut=-0.5,dt=1.,interp=True):
    td = np.where(np.max(d[:,:],axis=-1) > cut)
    ti = np.argmax(d[td,:] > th,axis=-1)[0]
    
    if interp:
        t = ti + (cut - d[td,ti]) / (d[td,ti] - d[td,ti-1]) # check ti isn't zero
    return dt*t[0], d.shape[0]

def extract_decay_bounds(d,th=-0.8,cut=-0.5,dt=1.):
    td = np.where(np.max(d[:,:],axis=-1) > cut)
    ti = np.argmax(d[td,:] > th,axis=-1)[0]

    return dt*ti, d.shape[0]

def get_new_decays(d,th=-0.8,cut=-0.5):
    """
    Compare ensembles of simulations at different resolutions or spectral cutoffs.
    Return the time-indices and trajectory indices for trajectories that decayed in one set of simulations but not the other.
    Requires the initial conditions to be "the same" in some way.

    Returns:
      - Decayed trajectories for each ensemble
      - Trajectory numbers that decayed in the "larger" simulation (smaller index)
      - Decay times of these trajectories
      - Trajectory numbers that decayed in "smaller" simulation (larger index)
      - Decay times of these trajectories
    """
    im = len(d)
    tc = np.min(np.array([dc.shape[0] for dc in d]))  # Only include common trajectories

    d_i = [ np.where(np.max(dc[:tc,:],axis=-1) > cut) for dc in d ]

    new_i= [ np.setdiff1d(d_i[i],d_i[i+1]) for i in range(im-1) ]
    t_i = [ np.argmax(d[i][new_i[i],:] > th, axis=-1) for i in range(im-1) ]

    new_ii = [ np.setdiff1d(d_i[i+1],d_i[i]) for i in range(im-1) ]
    t_ii = [ np.argmax(d[i+1][new_ii[i],:] > th, axis=-1) for i in range(im-1) ]
    return d_i, new_i, t_i, new_ii, t_ii
    
def quick_plot_len_scale(files,lv,t0=40.):
    dt = 1.5440808887540916
    t = [ decay_thresh(f,100) for f in files ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot( 2**(i-3)*(dt*np.sort(tc[0])-t0),np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax

# Make this more efficient with where statements
def decay_rate_from_survive(d,tmin,tmax):
    t = extract_decay_times(d)
    x = np.sort(t[0]); y = np.log(1.-np.linspace(0.,1.,t[1])[:x.size])
    imin = np.argmax(x>tmin); imax = np.argmin(x<tmax)
    return np.polyfit(x[imin:imax],y[imin:imax],1)

def decay_rate_from_decay(d,tmin,tmax):
    gamma = 1.
    return gamma

def t_off(th,cf,ct):
    return

if __name__=="__main__":
    base = 'runs/convergence/l1.2/len25/'
    files = [ base+'kcut0.25/n2048/bubble-count.dat',
              base+'kcut0.5/n1024/bubble-count.dat',
              base+'kcut1/n1024/bubble-count.dat',
              base+'kcut2/n1024/bubble-count.dat',
              base+'kcut4/n1024/bubble-count.dat'
              ]

    base = 'runs/rough/len12.5/'
    files_r = [ base+'l%s/n1024/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.55'] ]

    base = 'runs/convergence/spec_cut/len25/n2048/'; kv = ['2','3o8','4','3o16','8']
    files_k_1 = [ base+'kc%s/bubble-count.dat' % k for k in kv ]

    base = 'runs/convergence/spec_cut/len12.5/n1024/l1.2/'; kv = ['2','3o8','4','8','16']    
    files_k_2 = [ base+'kc%s/bubble-count.dat' % k for k in kv ]

    base = 'runs/convergence/spec_cut/len12.5/n2048/l1.2/'; kv = ['2','7o16','3o8','5o16','4']
    files_k_0 = [ base+'kc%s/bubble-count.dat' % k for k in kv ]

    base = 'runs/vary_phi0/l1.2/'
    files_p = [  base+'phi_%spi/bubble-count.dat' % s for s in ['0.3','0.35','0.4','0.45','0.5','0.55','0.6'] ]

    lab_tmp = ['1','2','3','4','5','6','7']

    base = 'runs/vary_const/phi0_0.4/'
    files_c_1 = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.3/'
    files_c_0 = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.5/'
    files_c_2 = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.4_kc3o8/'
    files_c_3 = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_lam/phi_0.3pi/'
    files_l = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.55','1.6'] ]

# phi0 variation
#    tmin = np.array([25.,35.,40.,50.,60.,])
#    tmax = np.array([35.,50.,60.,120.,200.,])

#    gammas = [ , -0.050195123 ]
#    y0 = [ , 1.14839629 ]
