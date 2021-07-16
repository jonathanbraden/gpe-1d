#include "macros.h"

!TO DO : Move fluctuation subroutines into a separate module
! Update thermal fluctuations, etc. to mimic the vacuum fluctuations (where I've fixed all of the horrible nonlocality)

program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use gaussianRandomField
  use eom
  use integrator
  use bubble_extraction!, only : count_bubbles, mean_cos, make_fourier_window, smooth_filter
!  use fluctuations
  
  implicit none
  real(dl), dimension(:,:), pointer :: fld
  real(dl), pointer :: time
  real(dl) :: dtout_, dt_

  integer, parameter :: inFile = 70, cpFile = 71
  integer :: i

  real(dl) :: alph, t_cross
  integer :: n_cross
  
  type SimParams
     real(dl) :: dx, dt, dtout
     integer :: nLat
  end type SimParams
  type(SimParams) :: sim
  
  fld(1:nLat,1:2) => yvec(1:2*nLat*nFld)
  time => yvec(2*nLat*nFld+1)

  alph = 8._dl; n_cross = 2
  
  call initialize_rand(93286123,12)
  call setup(nVar)

!  do i=1,200
!     call initialise_fields(fld,nLat/4)
!     call time_evolve(dx/alph,4*nlat*n_cross,128) ! Adjust this as needed!   
!     call time_evolve(0.4_dl/omega,10000,100)
!  enddo

  call initialise_fields(fld,nLat/2+1,0.25*twopi,(nLat/2+1)/2)
  call time_evolve(dx/alph,4*nlat*n_cross,128)

!  call forward_backward_evolution(0.4_dl/omega,10000,100)
  
contains

  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup
  !> before starting the program
  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
    call initialize_transform_1d(tPair,nLat)
  end subroutine setup

  subroutine initialise_fields(fld,kmax,phi,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax
    real(dl), intent(in), optional :: phi
    integer, intent(in), optional :: klat
    integer :: kc
    real(dl) :: phiL

    kc = nLat/2+1; if (present(klat)) kc = klat
    phiL = 0.5_dl*twopi; if (present(phi)) phiL = phi
    call initialise_mean_fields(fld)
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialize_vacuum_fluctuations(fld,len,m2eff,kmax,phiL,kc)
  end subroutine initialise_fields

  !>@brief
  !> Initialise Minkowski Gaussian vacuum approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Fix nonlocality of m2eff
  subroutine initialize_vacuum_fluctuations(fld,len,m2,kspec,phi0,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: len, m2
    integer, intent(in), optional :: kspec, klat
    real(dl), intent(in), optional :: phi0
    
    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spec, w2eff  ! remove w2eff here, it's unneeded
    real(dl), dimension(1:size(fld(:,1))) :: df
    integer :: i,km,kc
    real(dl) :: phiL, norm

    integer :: n, nn; real(dl) :: dk
    dk = twopi / len; n = size(fld(:,1)); nn = n/2+1
    
    km = size(spec); if (present(kspec)) km = kspec
    kc = size(spec); if (present(klat))  kc = klat
    
    phiL = twopi; if (present(phi0)) phiL = phi0

    ! The second 1/sqrt(2) is a bug since this factor isn't in my GRF sampler
    norm = (0.5_dl)**0.5 / phiL / sqrt(2._dl) / sqrt(len) ! second factor of 1/sqrt(2) is normalising the Box-Mueller, first one is from 1/sqrt(2\omega)

    do i=1,nn
       w2eff(i) = m2 + dk**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.25
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_vacuum_fluctuations

  !>@brief
  !> Initialise Minkowski Gaussian vacuum approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Fix nonlocality with len, m2eff, nlat, etc
  ! TO DO: Add an option to instead directly compare lattices of the same size with different spectral cuts
  subroutine initialize_thermal_fluctuations(fld,temp,kspec,phi0,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: temp
    integer, intent(in), optional :: kspec, klat
    real(dl), intent(in), optional :: phi0
    real(dl) :: df(1:nlat), spec(1:nLat/2+1), w2eff(1:nLat/2+1)
    integer :: i,km,kc, n
    real(dl) :: phiL, norm

    km = size(spec); if (present(kspec)) km = kspec
    kc = size(spec); if (present(klat))  kc = klat
    
    phiL = twopi; if (present(phi0)) phiL = phi0
    
    ! The 1/sqrt(2) is a bug since this factor isn't in my GRF sampler
    norm = (0.5_dl)**0.5 / phiL / sqrt(len) 

    do i=1,nLat/2+1
       w2eff(i) = m2eff + (twopi/len)**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.25 * sqrt( 1._dl/(exp(w2eff(2:km)**0.5/temp)-1._dl) + 0.5_dl)
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_thermal_fluctuations

  ! Add only thermal fluctuations
  !>@brief
  !> Initialise Minkowski Gaussian vacuum approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Fix nonlocality with len, m2eff, nlat, etc
  ! TO DO: Add an option to instead directly compare lattices of the same size with different spectral cuts
  subroutine initialize_only_thermal_fluctuations(fld,temp,kspec,phi0,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: temp
    integer, intent(in), optional :: kspec, klat
    real(dl), intent(in), optional :: phi0
    real(dl) :: df(1:nlat), spec(1:nLat/2+1), w2eff(1:nLat/2+1)
    integer :: i,km,kc, n
    real(dl) :: phiL, norm

    km = size(spec); if (present(kspec)) km = kspec
    kc = size(spec); if (present(klat))  kc = klat
    
    phiL = twopi; if (present(phi0)) phiL = phi0
    
    ! The 1/sqrt(2) is a bug since this factor isn't in my GRF sampler
    norm = (0.5_dl)**0.5 / phiL / sqrt(len) 

    do i=1,nLat/2+1
       w2eff(i) = m2eff + (twopi/len)**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.25 * sqrt( 1._dl/(exp(w2eff(2:km)**0.5/temp)-1._dl) )
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_only_thermal_fluctuations
  
  ! Add subroutine for high-temperature limit of thermal fluctuations
  !>@brief
  !> Initialise Minkowski Gaussian vacuum approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Fix nonlocality with len, m2eff, nlat, etc
  ! TO DO: Add an option to instead directly compare lattices of the same size with different spectral cuts
  subroutine initialize_high_T_fluctuations(fld,temp,kspec,phi0,klat)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: temp
    integer, intent(in), optional :: kspec, klat
    real(dl), intent(in), optional :: phi0
    real(dl) :: df(1:nlat), spec(1:nLat/2+1), w2eff(1:nLat/2+1)
    integer :: i,km,kc, n
    real(dl) :: phiL, norm

    km = size(spec); if (present(kspec)) km = kspec
    kc = size(spec); if (present(klat))  kc = klat
    
    phiL = twopi; if (present(phi0)) phiL = phi0
    
    ! The 1/sqrt(2) is a bug since this factor isn't in my GRF sampler
    norm = (0.5_dl)**0.5 / phiL / sqrt(len) 

    do i=1,nLat/2+1
       w2eff(i) = m2eff + (twopi/len)**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:km) = norm / w2eff(2:km)**0.5 * sqrt(temp)
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(:kc),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_high_T_fluctuations
  
  function light_cross_time(len) result(tmax)
    real(dl), intent(in) :: len
    real(dl) :: tmax
    tmax = 0.5_dl*len
  end function light_cross_time

  function convert_t_to_nstep(dt,dtout,tend) result(ns)
    real(dl), intent(in) :: dt,dtout, tend
    integer, dimension(1:2) :: ns

    ns(1) = int(tend/dt)
    ns(2) = int(dtout/dt)
  end function convert_t_to_nstep

  subroutine convert_tstep_to_int(dt,dtout,tend,ns,nout)
    real(dl), intent(in) :: dt, dtout, tend
    integer, intent(out) :: ns, nout

    ns = int(tend/dt)
    nout = int(dtout/dt)
  end subroutine convert_tstep_to_int

  !>@brief
  !> Initialise a mean field around which to sample fluctuations within the specified band.
  !> Returns the constrained part of the field, and the remaining fluctuations to sample
  subroutine constrained_fluctuations(imin,imax,ns)
    integer, intent(in) :: imin, imax, ns
    integer :: i
    real(dl), dimension(1:nlat/2+1) :: spec, spec_red
    
    ! Start by generating the "mean field"
    spec = 0._dl
    do i=1,size(spec)
       spec = 1._dl
    enddo
    spec_red = 0._dl; spec_red(imin:imax) = spec(imin:imax)
    ! Generate the field to sample around
    
    ! Now sample the non-constrained wavenumbers
    !!!!!!!!! Check the indexing in here to make sure I'm not accidentally resampling
    spec_red = 0._dl; spec_red(:imin) = spec(:imin); spec_red(imax:) = spec(imax:)
    do i=1,ns
       ! Sample the reduced spectrum
       ! time evolve
       ! compute any desired statistics
    enddo
  end subroutine constrained_fluctuations
  
  !>@brief
  !> Evolve a collection of ns field trajectories holding the long-wavelength part of the field fixed while varying the short wavelengths
  subroutine vary_high_k_modes(phi_l,ns)
    real(dl), dimension(:,:), intent(in) :: phi_l
    integer, intent(in) :: ns

    real(dl), dimension(1:nlat) :: df
    integer :: i
!    call initialise_fields(phi_l,nlat/8)
    do i=1,ns
       ! call generate_1dGRF(df)
       fld(:,1) = phi_l(:,1) + df
       ! call generate_1dGRF(df)
       fld(:,2) = phi_l(:,2) + df
       !call time_evolve()
    enddo
  end subroutine vary_high_k_modes

  !>@brief
  !> Evolve a collection of ns field trajectories holding the short-wavelength part of the field fixed while varying the long wavelengths
  subroutine vary_low_k_modes(phi_s,ns)
    real(dl), dimension(:), intent(in) :: phi_s
    integer, intent(in) :: ns

  end subroutine vary_low_k_modes
  
  subroutine forward_evolution(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no

    call initialise_fields(fld,nLat/8)
    call setup(nVar)
    call output_fields(fld)
    call time_evolve(dt,ns,no)
    call write_checkpoint(fld,time,cpFile)
  end subroutine forward_evolution
  
  subroutine forward_backward_evolution(dt,ns,no,amp)
    real(dl), intent(in) :: dt
    integer,intent(in) :: ns,no
    real(dl), intent(in), optional :: amp
    
!    call initialize_rand(72,18)
    call initialise_fields(fld,nLat/8)
    call setup(nVar)
    call output_fields(fld)
    call time_evolve(dt,ns,no)
    call write_checkpoint(fld,time,cpFile)

    ! now time reverse by flipping sign of time derivative
    fld(:,2) = -fld(:,2)
    if (present(amp)) call initialise_new_fluctuations(fld,amp)  

    call time_evolve(dt,ns,no)
  end subroutine forward_backward_evolution

  subroutine initialise_from_file(fName,n)
    character(*), intent(in) :: fName
    integer, intent(in) :: n
    integer :: i; real(dl) :: f,df
    integer :: fNum

    fNum = inFile
    open(unit=fNum,file=fName)
    do i=1,n
       read(fNum,*) f,df; fld(i,1) = f; fld(i,2) = df
    enddo
    close(fNum)
  end subroutine initialise_from_file

  !>@brief
  !> Reverse the time flow of a simulation by flipping the sign of phidot
  subroutine reverse_time(fld)
    real(dl), intent(inout), dimension(:,:) :: fld
    fld(:,2) = -fld(:,2)
  end subroutine reverse_time
  
  !>@brief
  !> Initialise the field to have mean value given by the false vacuum and no mean velocity 
  subroutine initialise_mean_fields(fld)
    real(dl), dimension(:,:), intent(out) :: fld
    fld(:,1) = 0.5_dl*twopi
    fld(:,2) = 0._dl
  end subroutine initialise_mean_fields

  subroutine initialise_fluc(fld,spec,kmax)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), dimension(:,:), intent(in) :: spec
    integer, intent(in), optional :: kmax

    real(dl) :: df(1:nLat); integer :: km

    km = size(spec(:,1)); if (present(kmax)) km = kmax
    call generate_1dGRF(df,spec(1:km,1),.false.)
    fld(:,1) = fld(:,1) + df
    call generate_1dGRF(df,spec(1:km,2),.false.)
    fld(:,2) = fld(:,2) + df
  end subroutine initialise_fluc
  
  !>@brief
  !> Initialise the field fluctuations
  subroutine initialise_fluctuations(fld,kmax,phi0)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in), optional :: kmax
    real(dl), intent(in), optional :: phi0
    
    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j, km
    real(dl) :: phiL
    
    km = size(spec); if (present(kmax)) km = kmax
    phiL = 0.5_dl; if (present(phi0)) phiL = phi0
    
    spec = 0._dl
    do i=2,nLat/2
       spec(i) = (1./phiL)*(0.5_dl)**0.5 / (sqrt(len*rho))
    enddo
    call generate_1dGRF(df,spec(1:km),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * m2eff**0.5
    call generate_1dGRF(df,spec(1:km),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialise_fluctuations

  !>@brief
  !> Initialise the field fluctuations
  subroutine initialise_new_fluctuations(fld,amp)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), intent(in) :: amp
    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j
    
    spec = 0._dl
    do i=2,nLat/2
       spec(i) = (0.5_dl)**0.5 / (sqrt(len*rho))
    enddo
    call generate_1dGRF(df,spec(1:128),.false.)
    fld(:,1) = fld(:,1) + amp*df(:)
    call generate_1dGRF(df,spec(1:128),.false.)
    fld(:,2) = fld(:,2) + amp*df(:)

  end subroutine initialise_new_fluctuations
  
  subroutine time_evolve(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no
    integer :: i,j, outsize, nums

    open(unit=60,file='bubble-count.dat')
    
    print*,"dx is ", dx, "dt is ",dt, "dx/dt is ",dx/dt
    if (dt > dx) print*,"Warning, violating Courant condition"
    
    outsize = ns/no; nums = ns/outsize
    print*,"dt out is ",dt*outsize
    dt_ = dt; dtout_ = dt*outsize
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_fields(fld)
       write(60,*) count_bubbles(fld(:,1),4), mean_cos(fld(:,1))
    enddo
    write(60,*)
  end subroutine time_evolve

#define SMOOTH 1
  subroutine output_fields(fld)
    real(dl), dimension(:,:), intent(in) :: fld
    logical :: o; integer :: i
    integer, parameter :: oFile = 99
    real(dl), dimension(1:nLat) :: gsq, gsq_fd
    
    real(dl) :: lambda
#ifdef SMOOTH
    real(dl), dimension(1:nLat) :: w_box, f_sm
    complex(dl), dimension(1:nLat/2+1) :: wk
#endif
    
    lambda = del*(2._dl/nu)**0.5
    
    inquire(file='fields.dat',opened=o)
    if (.not.o) then
       open(unit=oFile,file='fields.dat')
       write(oFile,*) "# Lattice Parameters"
       write(oFile,*) "# n = ",nLat," dx = ",dx
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_
    endif
    
    gsq_fd(1) = 0.5_dl*( (fld(nLat,1)-fld(1,1))**2+(fld(2,1)-fld(1,1))**2 )
    gsq_fd(nLat) = 0.5_dl*( (fld(nLat-1,1)-fld(nLat,1))**2+(fld(nLat,1)-fld(1,1))**2  )
    gsq_fd(2:nLat-1) = 0.5_dl*( (fld(1:nLat-2,1)-fld(2:nLat-1,1))**2+(fld(3:nLat,1)-fld(2:nlat-1,1))**2 )
    gsq_fd = gsq_fd / dx**2
#ifdef FOURIER
    tPair%realSpace(:) = fld(1:nLat,1)
    call gradsquared_1d_wtype(tPair,dk)
    gsq(:) = tPair%realSpace(:)
#else
    gsq(:) = 0._dl  ! tPair isn't created unless doing Fourier transforms
#endif

#ifdef SMOOTH
!    w_box = 0._dl
!    w_box(1:51) = 1._dl; w_box(nLat-50+1:nLat)=1._dl
    call box_filter_r(nLat,50,w_box)
    call make_fourier_window(w_box,wk,tPair)
    call smooth_filter(f_sm,fld(:,1),wk,tPair)
#endif
    
    ! Fix this if I change array orderings
    do i=1,nLat
       write(oFile,*) fld(i,:), gsq(i), 4._dl*nu*(-cos(fld(i,1)) + 0.5_dl*lambda**2*sin(fld(i,1))**2 - 1._dl), gsq_fd(i), 2._dl*nu*(-1.+lambda**2)*(fld(i,1)-0.5_dl*twopi)**2, f_sm(i)
    enddo
    write(oFile,*)
    
  end subroutine output_fields

  !>@brief
  !> Write a checkpoint file with all information for restarting the simulation.
  subroutine write_checkpoint(fld,tcur,fn)
    real(dl), intent(in) :: fld(:,:), tcur
    integer, intent(in) :: fn
    integer :: i
    open(unit=fn,file='flds.chk')
    write(fn,*) nLat, dx, tcur
    do i=1,nLat
       write(fn,*) fld(i,:)
    enddo
    close(fn)
  end subroutine write_checkpoint

  !>@brief
  !> Read in a previously produced checkpoint file to initialise the simulation.
  subroutine read_checkpoint(fld,tcur,nLat,fName)
    real(dl), intent(out) :: fld(:,:), tcur
    integer, intent(out) :: nLat
    character(*), intent(in) :: fName

    integer, parameter :: fn = 73
    integer :: i
    integer :: n
    real(dl) :: dx, tc
    
    open(unit=fn,file=fName)  ! Fix nonlocality
    read(fn,*) n, dx, tc
    print*,"Reading checkpoint file: N = ",n," dx = ",dx," t = ",tc
    ! Add a check that n matches the parameter used for the sim
    do i=1,n
       read(fn,*) fld(i,:)
    enddo
    close(fn)
  end subroutine read_checkpoint
  
end program Gross_Pitaevskii_1d
