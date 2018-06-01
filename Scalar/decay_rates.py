#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def get_trajectories(files,nt,ax=-1):
    d = []
    for f in files:
        a = np.loadtxt(f)
        d.append( a[:a.shape[0]-a.shape[0]% nt,ax].reshape((-1,nt)) )
    return d

def quick_plot_traj(d,lv,**kwargs):
    t = [ extract_decay_times(dc,th=t,cut=tc,interp=interp) for dc in d ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm decay}$')
    ax.set_xlabel(r'$\bar{t}$')
    ax.set_ylim(-0.1,1.1)

    return fig,ax

def quick_plot_traj_survive(d,lv,dt=1.,**kwargs):
    t = [ extract_decay_times(dc,**kwargs) for dc in d ]

    fig, ax = plt.subplots()
    for i,tc in enumerate(t):
        ax.plot(dt*np.sort(tc[0]),1.-np.linspace(0.,1.,tc[1])[:tc[0].size],label=lv[i])
    plt.legend(loc='lower right',ncol=2)
    ax.set_ylabel(r'$P_{\rm survive}$')
    ax.set_xlabel(r'$\bar{t}$')
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

def extract_decay_times(d,th=-0.8,cut=-0.5,interp=False,dt=1.,**kwargs):
    td = np.where(np.max(d[:,:],axis=-1) > cut)
    ti = np.argmax(d[td,:] > th,axis=-1)[0]
    
    # This needs to be fixed to work when ti = 0, or fucky slope
    # Why does this thing fuck up so bad?  Because if I'm using early times, might not have an increasing slope, which will fuck everything up ...
    if interp:
        t = ti + np.sign(ti)*(th - d[td,ti]) / (d[td,ti] - d[td,ti-1]) 
        t = t[0]
    else:
        t = ti
    return dt*t, d.shape[0]

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

def plot_compare(d_1,d_2,direct=True,**kwargs):
    if (len(d_1) != len(d_2)):
        print("Error, incompatible lengths")
        return

    fig,ax = plt.subplots()
    for i in range(len(d_1)):
        im1 = d_1[i].shape[0]; im2 = d_2[i].shape[0]
        if (direct):
            im = min(im1,im2); im1=im; im2=im
        t = extract_decay_times(d_1[i][:im1],**kwargs)
        ax.plot(np.sort(t[0]),1.-np.linspace(0.,1.,t[1])[:t[0].size])
        t = extract_decay_times(d_2[i][:im2],**kwargs)
        ax.plot(np.sort(t[0]),1.-np.linspace(0.,1.,t[1])[:t[0].size],'--')
    return fig,ax

def plot_thresh_scan(d,tv,**kwargs):
    fig,ax = plt.subplots()
    for th in tv:
        t = extract_decay_times(d,th=th,**kwargs)
        ax.plot(np.sort(t[0]),1.-np.linspace(0.,1.,t[1])[:t[0].size])
    return fig,ax

def plot_compare_interp(d,tv,**kwargs):
    fig,ax = plt.subplots()
    t = extract_decay_times(d,th=tv,interp=False,**kwargs)
    ax.plot(np.sort(t[0]),1.-np.linspace(0.,1.,t[1])[:t[0].size])
    t = extract_decay_times(d,th=tv,interp=True,**kwargs)
    ax.plot(np.sort(t[0]),1.-np.linspace(0.,1.,t[1])[:t[0].size],'--')
    return fig, ax

def plot_lin_fits():
    return

def survival_plot(t,ax,style='-'):
    ax.plot(np.sort(t[0]),1.-np.linspace(0.,1.,t[1])[:t[0].size],style)
    return

def get_tbounds_from_frac(d,fmin,fmax,**kwargs):
    t = extract_decay_times(d,**kwargs)
    y = 1.-np.linspace(0.,1.,t[1])
    x = np.sort(t[0])
    imin = np.argmax(y<fmax); imax = np.argmax(y<fmin)
    tmin = x[imin]; tmax = x[imax]
    return tmin,tmax

def decay_rates_vary_frac(d,fmin,fmax,**kwargs):
    t = [ extract_decay_times(dc,**kwargs) for dc in d ]
    gam = np.empty((len(d),len(fmin)))
    for j,(f0,f1) in enumerate(zip(fmin,fmax)):
        tb = [ get_tbounds_from_frac(dc,f0,f1,**kwargs) for dc in d ]
        gam[:,j] = [ decay_rate_from_survive(d[i],tb[i][0],tb[i][1],**kwargs)[0][0] for i in range(len(d)) ]
    return -gam

def decay_rates_vary_thresh(d,tv,fmin,fmax,**kwargs):
    gam = np.empty((tv.size,len(d)))
    for j,th in enumerate(tv):
        t = [ extract_decay_times(dc,th=th,cut=th,**kwargs) for dc in d ] # fix cut
        tb = [ get_tbounds_from_frac(dc,fmin,fmax,th=th,cut=th,**kwargs) for dc in d ]
        gam[j,:] = [ decay_rate_from_survive(d[i],tb[i][0],tb[i][1],th=th,cut=th,**kwargs)[0][0] for i in range(len(d)) ]
    return -gam

def scan_decay_rates(d,tmin,tmax,**kwargs):
    p = np.array([ decay_rate_from_survive(d[i],tmin[i],tmax[i],**kwargs)[0] for i in range(len(d)) ])
    x = [ decay_rate_from_survive(d[i],tmin[i],tmax[i],t)[1] for i in range(len(d)) ]
    y = [ decay_rate_from_survive(d[i],tmin[i],tmax[i],t)[2] for i in range(len(d)) ]
    return p,x,y

# Make this more efficient with where statements
def decay_rate_from_survive(d,tmin,tmax,**kwargs):
    t = extract_decay_times(d,**kwargs)
    x = np.sort(t[0]); y = np.log(1.-np.linspace(0.,1.,t[1])[:x.size])
    imin = np.argmax(x>tmin); imax = np.argmin(x<tmax)
    return np.polyfit(x[imin:imax],y[imin:imax],1), x, y

def decay_rate_from_frac(d,fmin,fmax,**kwargs):
    t = extract_decay_times(d,**kwargs); y = 1.-np.linspace(0.,1.,t[1])
    x = np.sort(t[0])
    imin = np.argmax(y<fmax); imax = np.argmax(y<fmin)
    y = np.log(y[:x.size])
    return np.polyfit(x[imin:imax],y[imin:imax],1), x, y
    
def get_single_bubbles(d):
    """
    Set a constraint on the allowed slope to try and only obtain single bubble nucleations
    """
    return

def get_fv_rms(d):
    """
    Get the extraction threshold based on an empirical determination of the RMS
    """
    return np.array([np.std(dc[:,0]) for dc in d])

def instanton_rate(act,phi):
    return np.exp(-2.*np.pi*phi**2*act)*2.*np.pi*act*phi**2

def instanton_rate_0th(act,phi):
    """
    Decay rate ignoring scaling behaviour of determinant.  Normalised to unity when phi=1
    """
    return np.exp(-2.*np.pi*phi**2*act)*2.*np.pi*act

def create_binary(fName,files,n,dt,dx):
    d = get_trajectories(files,n)
    np.savez(fName,dt=np.array([dt]),dx=np.array([dx]))
    return

def merge_directories(fList,n,baseFile='bubble-count.dat'):
    d = get_trajectories(fList,n)
    return np.vstack(d)  # check this versus hstack, etc

def merge_L50_runs(dt,dx,kv='2'):
    base = 'runs/vary_phi0/l1.2/len50/kcut'+kv+'/'
    d = get_trajectories([base+'phi_%spi/bubble-count.dat' %s for s in ['0.25','0.3','0.35']],64)
    d2 = get_trajectories([base+'phi_%spi/bubble-count.dat' %s for s in ['0.4','0.45']],128)
    d65 = merge_directories([base+'phi_0.65pi%s/bubble-count.dat' %s for s in ['','_run2','_run3','_run4']],512)
    d60 = merge_directories([base+'phi_0.6pi%s/bubble-count.dat' %s for s in ['','_run2','_run3','_run4']],512)
    d55 = merge_directories([base+'phi_0.55pi%s/bubble-count.dat' %s for s in ['','_run2','_run3','_run4']],512)
    d50 = merge_directories([base+'phi_0.5pi%s/bubble-count.dat' %s for s in ['','_run2']],256)
    np.savez('cos-traj-vary-phi0-len50-kc'+kv+'-n4096.npz',dt=[dt],dx=[dx],p25=d[0],p30=d[1],p35=d[2],p40=d2[0],p45=d2[1],p50=d50,p55=d55,p60=d60,p65=d65)

if __name__=="__main__":
    dt = 1.5440808887540916; dx = 0.19301011109426144

    base_f = 'bubble-count.dat'
    lab_tmp = ['1','2','3','4','5','6','7','8','9']

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
    files_k_1 = [ base+'kc%s/'%k+base_f for k in kv ]

    base = 'runs/convergence/spec_cut/len12.5/n1024/l1.2/'; kv = ['2','3o8','4','8','16']    
    files_k_2 = [ base+'kc%s/'%k+base_f  for k in kv ]

    base = 'runs/convergence/spec_cut/len12.5/n2048/l1.2/'; kv = ['2','7o16','3o8','5o16','4']
    files_k_0 = [ base+'kc%s/bubble-count.dat' % k for k in kv ]


#    base = 'runs/vary_phi0/l1.2/len25/kcut4_old/'   # 512
#    files_p_0_a_25 = [  base+'phi_%spi/bubble-count.dat' % s for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65'] ]
    base = 'runs/vary_phi0/l1.2/len25/kcut4/'
    files_p_0_a_25 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]
    base = 'runs/vary_phi0/l1.2/len25/kcut3o8/' # 256
    files_p_0_b_25 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]
    base = 'runs/vary_phi0/l1.2/len25/kcut2/'   # 256
    files_p_0_c_25 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]

    base = 'runs/vary_phi0/l1.2/len12.5/kcut4/'
    files_p_0_a_12 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]  # steps are [128,128,128,128,128,128,128,256]
    base = 'runs/vary_phi0/l1.2/len12.5/kcut3o8/'
    files_p_0_b_12 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]  # steps are [128,128,128,128,128,128,128,256]
    base = 'runs/vary_phi0/l1.2/len12.5/kcut2/'
    files_p_0_c_12 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]  # steps are [128,128,128,128,128,128,128,256]
    base = 'runs/vary_phi0/l1.2/len12.5/kcut4_n2048/'
    files_p_0_cx_12 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]  # steps are [128,128,128,128,128,128,128,256]
    base = 'runs/vary_phi0/l1.2/len12.5/kcut2_n2048/'
    files_p_0_dx_12 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']]

    base = 'runs/vary_phi0/l1.2/len50/kcut2_old/'
    files_p_0_c_50 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6'] ] # Up to 0.35 is 128, after is 512
    #base = 'runs/vary_phi0/l1.2/len50/kcut2/'
    #files_p_0_c_50 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65'] ] # steps are [64,64,64,128,128,256,512,512,512]
    base = 'runs/vary_phi0/l1.2/len50/kcut3o8/'
    files_p_0_b_50 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65'] ] # steps are [64,64,64,128,128,256,512,512,512
    base = 'runs/vary_phi0/l1.2/len50/kcut4/'
    files_p_0_a_50 = [ base+'phi_%spi/'%s+base_f for s in ['0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6','0.65']] # steps are [64,64,64,128,128,256,512,512,512

    base = 'runs/vary_phi0/l1.3/'
    files_p_1_a = [ base+'phi_%spi/'%s+base_f for s in ['0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5','0.55','0.6']] # 512 steps per sim

    base = 'runs/vary_const/phi0_0.4/'
    files_c_1_a = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.4_kc3o8/'
    files_c_1_b = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.4_kc2/'
    files_c_1_c = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_const/phi0_0.3/'
    files_c_0_a = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.3_kc3o8/'
    files_c_0_b = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]
    base = 'runs/vary_const/phi0_0.3_kc2/'
    files_c_0_c = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_const/phi0_0.5/'
    files_c_2 = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_lam/phi_0.3pi/'
    files_l_0_a = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_lam/phi_0.3pi/kcut3o8/'
    files_l_0_b = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_lam/phi_0.3pi/kcut2/'
    files_l_0_c = [ base+'l%s/bubble-count.dat' % s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

    base = 'runs/vary_lam/phi_0.4pi/'
    files_l_1 = [ base+'l%s/bubble-count.dat' %s for s in ['1.2','1.3','1.4','1.5','1.6'] ]

# phi0 variation
    phi = np.sqrt(0.5)*np.pi*np.array([0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6])

    # These are for k_nyq/2 cut
    tmin = np.array([25.,35.,40.,50.,60.,100.,100.])
    tmax = np.array([35.,50.,80.,120.,200.,400.,400.])

    act = np.array([0.7299,1.0897,1.45,1.8155,2.0002,2.1869])

    # These are taken off the k_nyq cut, but will be used everywhere to start with, subject to fixing
    fmin = np.array([0.05,0.05,0.05,0.05,0.05,0.5,0.8,0.94])
    fmax = np.array([0.2,0.2,0.2,0.3,0.4,0.7,0.9,0.97])
    thresh=-0.5; cut=-0.5
    fD = ['p25','p30','p35','p40','p45','p50','p55','p60']
### L=25, k_cut=k_nyq

    dFile = 'runs/vary_phi0/l1.2/len25/kcut2/cos-traj-vary-phi0-len25-kc2-n2048.npz'
    d_k2_25 = np.load(dFile); a_c = [ d_k2_25[s] for s in fD ]
    gam_k2_l25 = -np.array( [ decay_rate_from_frac(a_c[i],fmin[i],fmax[i],th=thresh,cut=cut,interp=True)[0][0] for i in range(len(a_c)) ])

### L=25, k_cut=3*k_nyq/4    
    dFile = 'runs/vary_phi0/l1.2/len25/kcut3o8/cos-traj-vary-phi0-len25-kc3o8-n2048.npz'
    d_k3o8_25 = np.load(dFile); a_b = [ d_k3o8_25[s] for s in fD ]
    gam_k3o8_l25 = -np.array( [ decay_rate_from_frac(a_b[i],fmin[i],fmax[i],th=thresh,cut=cut,interp=True)[0][0] for i in range(len(a_b)) ])

### L=25, k_cut=k_nyq/2
    dFile = 'runs/vary_phi0/l1.2/len25/kcut4/cos-traj-vary-phi0-len25-kc4-n2048.npz'
    d_k4_25 = np.load(dFile); a_a = [ d_k4_25[s] for s in fD ]
    gam_k4_l25 = -np.array( [ decay_rate_from_frac(a_a[i],fmin[i],fmax[i],th=thresh,cut=cut,interp=True)[0][0] for i in range(len(a_a))])

#### L = 12.5 ##########
    fmin = np.array([0.05,0.05,0.05,0.2,0.53,0.84,0.96])
    fmax = np.array([0.2,0.2,0.2,0.3,0.65,0.92,0.98])
    thresh=-0.5; cut=-0.5
    fD = ['p25','p30','p35','p40','p45','p50','p55']
### L=12.5, k_cut=2
    dFile = 'runs/vary_phi0/l1.2/len12.5/kcut2/cos-traj-vary-phi0-len12.5-kc2-n1024.npz'
    d_k2_12 = np.load(dFile); a_a = [ d_k2_12[s] for s in fD ]
    gam_k2_l12 = -np.array( [ decay_rate_from_frac(a_a[i],fmin[i],fmax[i],th=thresh,cut=cut,interp=True)[0][0] for i in range(len(a_a))])
### Add the other spectral cuts here when they finish running

    fD = ['p25','p30','p35','p40','p45','p50','p55','p60','p65']
    fmin = np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.5,0.5,0.96])
    fmax = np.array([0.2,0.2,0.2,0.2,0.2,0.2,0.7,0.8,0.98])
### L = 50, k_cut=2
    dFile = 'runs/vary_phi0/l1.2/len50/kcut2_old/cos-traj-vary-phi0-len50-kc2-n4096.npz'
    d = np.load(dFile); a = [ d[s] for s in fD ]
    gam_k2_l50 = -np.array( [ decay_rate_from_frac(a[i],fmin[i],fmax[i],th=thresh,cut=cut,interp=True)[0][0] for i in range(len(a))])
### Add the other spectral cuts and fix the above one


#    p_0,x_0,y_0 = scan_decay_rates(a,tmin,tmax,-0.5)
    

#    d_p = get_trajectories(files_p,512)
#    p_0,x_0,y_0 = scan_decay_rates(d_p,tmin,tmax,-0.5)
#    p_1,x_1,y_1 = scan_decay_rates(d_p,tmin,tmax,-0.6)

#    d_l = get_trajectories(files_l,512)
#    p_l,x_l,y_l = scan_decay_rates(d_l,25.*np.ones(len(d_l)),40.*np.ones(len(d_l)),-0.5)
