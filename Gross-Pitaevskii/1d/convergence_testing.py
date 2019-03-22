#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# Define some notation
t_lab = r'$\bar{t}$'
norm_lab = r'$||_{\rm max}$'

files = [ 'paper-data/fields-n256-w250-l1.4-dt%s.dat' % s for s in ['1','2','4','8'] ]

files_floquet_varydt_1024 = [ 'paper-data-new/convergence/floquet-dt/fields-floquet-n1024-%swosc.dat' %s for s in [ '8', '16', '32', '64', '128', '256' ] ]
files_floquet_varydt_2048 = [ 'paper-data-new/convergence/floquet-dt/fields-floquet-n2048-%swosc.dat' %s for s in [ '32', '64', '128', '256' ] ]
files_floquet_varydt_290 = [ 'paper-data-new/convergence/floquet-dt/fields-floquet-n290-%swosc.dat' %s for s in [ '2','4','8','16','32','64','128' ] ]
files_floquet_varydt_295 = [ 'paper-data-new/convergence/floquet-dt/fields-floquet-n295-%swosc.dat' %s for s in [ '2','4','8','16','32','64','128' ] ]
files_floquet_varydt_312 = [ 'paper-data-new/convergence/floquet-dt/fields-floquet-n312-%swosc.dat' %s for s in [ '2','4','8','16','32','64','128'] ]

files_floquet_varydx_295_w128 = [ 'paper-data-new/convergence/floquet-dx/fields-floquet-n%s-w128.dat' %s for s in [ '295', '590', '1180', '2360' ] ]
files_floquet_varydx_295_w256 = [ 'paper-data-new/convergence/floquet-dx/fields-floquet-n%s-w256.dat' %s for s in [ '295', '590', '1180', '2360', '4720' ] ]

files_floquet_varydt_discrete_512 = [ 'paper-data-new/convergence/floquet-dt-discrete/fields-floquet-n512-w%s.dat' %s for s in [ '4','8','16','32','64','128', '256' ] ]

files_no_floquet_varydx_512_w256 = [ 'paper-data-new/convergence/no-floquet-dx/fields-floquet-n{0:s}-wosc{1:s}.dat'.format(n,w) for n,w in zip([ '512','1024','2048','4096'],['256','256','256','256']) ]

files_no_floquet_varydx_295 = [ 'paper-data-new/convergence/no-floquet-dx/fields-floquet-n{0:s}-w{1:s}.dat'.format(n,w)
                                for n,w in zip([ '295', '590', '1180', '2360', '4720' ],
                                ['128','128','256','256','256']) ]

files_no_floquet_varydx_dtlong = [ 'paper-data-new/convergence/no-floquet-dx/fields-floquet-n{0:s}-w{1:s}.dat'.format(n,w) for n,w in zip([ '295','590','1180','2360','4720'],['128','128','128','128','256']) ]

files_no_floquet_varydt_512 = [ 'paper-data-new/convergence/no-floquet-dt/fields-floquet-n512-wosc%s.dat' %s for s in [ '4','8','16','32','64' ] ]

def collectDataFixedN(files,n):
    return np.array( [ np.genfromtxt(f).reshape((-1,n,4)) for f in files ] )

# nMin is number of grid points on smallest grid.  Assumed to increase by factor of 2 each time
def collectDataFixedT(files,nMin):
    return np.array([ np.genfromtxt(f).reshape((-1,nMin*2**ns,4))[:,::2**ns,:] for ns,f in enumerate(files) ] )

def plot_t_convergence_axis(a,d,dt=1.):
    tv = dt*np.arange(d.shape[1])
    for i in range(d.shape[0]-1):
        a.plot(tv,np.max(np.max(np.abs(d[i+1]-d[i]),axis=-1),axis=-1))
    a.set_yscale('log')
    a.set_xlim((tv[0],tv[-1]))
    return

def plot_x_convergence_axis(a,d,dt=1.):
    tv = dt*np.arange(d.shape[1])
    for i in range(d.shape[0]-1):
        a.plot(tv,np.max(np.max(np.abs(d[i+1]-d[i]),axis=-1),axis=-1))
    a.set_yscale('log')
    a.set_xlim((tv[0],tv[-1]))
    return
        
def plot_t_convergence(d,dt=1.):
    tv = dt*np.arange(d.shape[1])
    f,a = plt.subplots()
    for i in range(d.shape[0]-1):
        a.plot(tv,np.max(np.abs(d[i+1,:,:,0]-d[i,:,:,0]),axis=-1))
    a.set_yscale('log')
    a.set_xlabel(t_lab); a.set_ylabel(norm_lab)
    a.set_xlim((tv[0],tv[-1]))
    return f,a

def plot_t_convergence_rho_tot(d):
    f,a = plt.subplots()
    for i in range(d.shape[0]-1):
        a.plot(np.max(np.abs(np.sum(d[i+1]**2,axis=-1)-np.sum(d[i]**2,axis=-1)),axis=-1))
    a.set_yscale('log')
    a.set_xlabel(t_lab); a.set_ylabel(norm_lab)
    return f,a

def plot_x_convergence(d,dt=1.):
    f,a = plt.subplots()
    tv = dt*np.arange(d.shape[1])
    for i in range(d.shape[0]-1):
        a.plot(tv,np.max(np.abs(d[i+1,:,:,0]-d[i,:,:,0]),axis=-1))
    a.set_yscale('log')
    a.set_xlabel(t_lab); a.set_ylabel(norm_lab)
    return f,a

def n512_dt_convergence_plot():
    d = collectDataFixedN(files_varyt_n512,512)
    f,a = plot_t_convergence(d)
    f.show()
    f.savefig('gpe-convergence-dt-n512.pdf')
    return

def n256_dt_convergence_plot():
    d = collectDataFixedN(files_varyt_n256,256)
    f,a = plot_t_convergence(d)
    f.show()
    f.savefig('gpe-convergence-dt-n256.pdf')
    return

def dx_convergence_plot():
    d = collectDataFixedT(files_varyx_fix_w32,128)
    f,a = plot_x_convergence(d)
    return f,a

def dx_no_floquet_plot():
    d = collectDataFixedT(files_no_floquet_varydx_w256,512)
    f,a = plot_x_convergence(d)
    return f,a

def dt_no_floquet_plot():
    d = collectDataFixedN(files_no_floquet_varydt_512,512)
    f,a = plot_t_convergence(d)
    return f,a
    
def main():
    return

if __name__=="__main__":
    pass
