#include "macros.h"

program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use gaussianRandomField
  use eom
  use integrator
  use bubble_extraction, only : count_bubbles, mean_cos

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

  alph = 4._dl; n_cross = 2

  call initialize_rand(87,18)
  call setup(nVar)
  do i=1,200
     call initialise_fields(fld,nLat/8)
     call time_evolve(dx/alph,4*nlat*n_cross),128) ! Adjust this as needed
!     call time_evolve(0.4_dl/omega,10000,100)
  enddo
 
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

  subroutine initialise_fields(fld,kmax)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: kmax
    integer :: i; real(dl) :: dt, theta
    integer :: km
    
    call initialise_mean_fields(fld)
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here
    !    call initialise_fluctuations(fld)
    call initialize_vacuum_fluctuations(fld,kmax)
  end subroutine initialise_fields
  
  subroutine convert_tstep_to_int(dt,dtout,tend,ns,nout)
    real(dl), intent(in) :: dt, dtout, tend
    integer, intent(out) :: ns, nout

    ns = int(tend/dt)
    nout = int(dtout/dt)
  end subroutine convert_tstep_to_int
  
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

    print*,"Mean field is ", sum(fld(:,1))/nLat - 0.5_dl*twopi
    print*,"Mean dfld is ", sum(fld(:,2))/nLat
  end subroutine initialise_fluctuations

  !>@brief
  !> Initialise Minkowski Gaussian vacuum approximation for fluctuations.
  !> Spectra in this subroutine are truncated for direct comparison of 
  !> fluctuations generated between lattices of varying size.
  !
  ! TO DO: Fix nonlocality with len, m2eff, nlat, etc
  ! TO DO: Add an option to instead directly compare lattices of the same size with different spectral cuts
  subroutine initialize_vacuum_fluctuations(fld,kmax,phi0)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in), optional :: kmax
    real(dl), intent(in), optional :: phi0
    real(dl) :: df(1:nlat), spec(1:nLat/2+1), w2eff(1:nLat/2+1)
    integer :: i,km, n
    real(dl) :: phiL, norm

    km = size(spec); if (present(kmax)) km = kmax
    phiL = twopi; if (present(phi0)) phiL = phi0
    
    norm = (0.5_dl)**0.5 / phiL / sqrt(2._dl) / sqrt(len) ! second factor of 1/sqrt(2) is normalising the Box-Mueller, first one is from 1/sqrt(2\omega)
    ! Fix the nonlocality of the length

    do i=1,nLat/2+1
       w2eff(i) = m2eff + (twopi/len)**2*(i-1)**2
    enddo
    spec = 0._dl
    spec(2:) = norm / w2eff(2:)**0.25
    call generate_1dGRF(df,spec(1:km),.false.)
    fld(:,1) = fld(:,1) + df(:)

    spec = spec * w2eff**0.5
    call generate_1dGRF(df,spec(1:km),.false.)
    fld(:,2) = fld(:,2) + df(:)
  end subroutine initialize_vacuum_fluctuations
  
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
       write(60,*) count_bubbles(fld(:,1)), mean_cos(fld(:,1))
    enddo
    write(60,*)
  end subroutine time_evolve

  subroutine output_fields(fld)
    real(dl), dimension(:,:), intent(in) :: fld
    logical :: o; integer :: i
    integer, parameter :: oFile = 99
    real(dl), dimension(1:nLat) :: gsq, gsq_fd

    real(dl) :: lambda
    lambda = del*(2._dl/nu)**0.5

    if (.true.) return
    
    inquire(file='fields.dat',opened=o)
    if (.not.o) then
       open(unit=oFile,file='fields.dat')
       write(oFile,*) "# Lattice Parameters"
       write(oFile,*) "# n = ",nLat," dx = ",dx
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_
!       write(oFile,*) "# Model Parameters"
!       write(oFile,*) "# nu = ",nu," g = ",gs, " w = ",omega
!       write(oFile,*) "# rho = ", rho
!       write(oFile,*) "# delta = ", del
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
    ! Fix this if I change array orderings
    do i=1,nLat
       write(oFile,*) fld(i,:), gsq(i), 4._dl*nu*(-cos(fld(i,1)) + 0.5_dl*lambda**2*sin(fld(i,1))**2 - 1._dl), gsq_fd(i), 2._dl*nu*(-1.+lambda**2)*(fld(i,1)-0.5_dl*twopi)**2
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
