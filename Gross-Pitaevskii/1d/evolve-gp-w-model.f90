program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use utils, only : newunit
  use gaussianRandomField
  use eom
  use TimeStepping
  use homogeneousField
  use fluctuations
  use output
  use integrator

  implicit none

  type Model
     !type(ModelParams) :: model_params
     type(SpecParams) :: spec_params
     type(TimeStepper) :: time_params
  end type Model
  
  real(dl), dimension(:,:,:), pointer :: fld
  real(dl), pointer :: tcur
  real(dl) :: dtout_, dt_  ! Remove this ugliness

  integer :: n
  integer :: ds, w_samp

  ! Fixing parameters for now
  real(dl), parameter :: omega_s = 50._dl, rho_s = 1000._dl, del_s = 1._dl+0.4_dl
  real(dl), parameter :: len_s = 100._dl

  real(dl) :: phi0, kfloq
  integer :: wn

  integer :: i, nSamp

  ! Move params somewhere else
  real(dl) :: om_alex, rho_alex, nu_alex, del_alex, phi_init, len_alex

  type(SpecParams) :: fluc_params
  type(TimeStepper) :: time_stepper
  
  n=512
  ds = 1; w_samp=16

  ! Old stuff for false vacuum
  !call initialise_model(2, del_s, 2.e-3, 1._dl, 0.5_dl, 0.5_dl, omega_s, rho_s)
  !call setup_w_model( 2, n, len_s/(2.*(nu)**0.5)/dble(n) )
  !fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  !tcur => yvec(2*nLat*nFld+1)
  !call initialise_false_vacuum(fld)
  !call time_evolve((1._dl/dble(w_samp))*twopi/omega,100._dl,dble(ds)/dble(w_samp)*twopi/omega)

  
  ! Turn these into parameters for readability
  len_alex = 454.65 
  n = 512  !1024
  w_samp = 16
  rho_alex = 1000. !43.99 !1000. !43.99 !43.99 
  nu_alex = 0.01  !0.01
  del_alex = 1.2_dl !1.2_dl
  om_alex = 64./2./sqrt(nu_alex)
  
  ! This is for sampling the ICs, which seems to be working
  !call initialise_model_symmetric(2, nu_alex, del_alex, om_alex, rho_alex)
  !call setup_w_model( 2, n, len_alex/n )
  !call sample_ics(1000, n, 2, n/4, 100.)
  
  phi_init = 0._dl ! 0.01*twopi
  call initialise_model_symmetric(2, nu_alex, del_alex, om_alex, rho_alex) ! fix rho
  call setup_w_model( 2, n, len_alex/n )
  fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  tcur => yvec(2*nLat*nFld+1)

  fluc_params%rho = rho_alex
  fluc_params%len = len_alex
  fluc_params%type = 'KG'
  fluc_params%modes = (/.true.,.true./)
  fluc_params%nCut = n/2
  fluc_params%nu = nu_alex
  fluc_params%m2eff = 4._dl*nu_alex*(del_alex**2-1._dl)
  
  call sample_ics(1000,fluc_params, n, 2)
  
  ! Initialise the fluctuations
  fld = 0._dl
  call initialise_fluctuations(fld, fluc_params)
  !call initialise_fluctuations_white(fld, rho_alex, nCut=n/2)
  !call initialise_fluctuations_bogoliubov(fld, rho_alex, nCut=n/2)
  !call initialise_fluctuations_kg_long(fld, rho, 4*nu*(del_alex**2-1._dl), n/2, (/.true.,.false./))
  !call initialise_fluctuations_white(fld, fluc_params)
  ! Add displaced mean field
  fld(:,1,1) = fld(:,1,1) + 1.
  fld(:,1,2) = fld(:,1,2) - 1.
  !fld(:,1,1) = cos(0.5*phi_init) + fld(:,1,1); fld(:,2,1) = -sin(0.5*phi_init) + fld(:,2,1)
  !fld(:,1,2) = -cos(0.5*phi_init) + fld(:,1,2); fld(:,2,2) = -sin(0.5*phi_init) + fld(:,2,2)

  ! Fix this up, then change the time-stepper call
  time_stepper%dt = (twopi/om_alex)/dble(w_samp)
  time_stepper%out_size = 8  ! Fix this
  time_stepper%tcur = 0._dl
  time_stepper%n_out_steps = 2.*45.456 / ((twopi/om_alex)/dble(2.))
  !call time_evolve( (twopi/om_alex)/dble(w_samp), 2.*45.465, (twopi/om_alex)/dble(2.) )
  !call time_evolve( (twopi/om_alex)/dble(w_samp), 

  !call time_evolve( (twopi/om_alex)/dble(w_samp), 45.465, (twopi/om_alex)/dble(w_samp)*8 )
  !call time_evolve( (twopi/om_alex)/dble(w_samp), sqrt(0.44)*twopi*2.*sqrt(nu_alex)*5, (twopi/om_alex)/dble(w_samp)*8 )
  ! For testing homogeneous dynamics
  !call time_evolve( (twopi/om_alex)/dble(w_samp), twopi/sqrt(del_alex**2-1.)/(2.*sqrt(nu))*8., twopi/sqrt(del_alex**2-1.)/(2.*sqrt(nu))/32. )

  
#ifdef PREHEAT_SINGLE
  call initialise_model_symmetric(2, 1.e-2, 0._dl, 0._dl, 1000.)
  call setup_w_model(2, n, len_s/(2.*(nu)**0.5)/dble(n))
  fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  tcur => yvec(2*nLat*nFld+1)
  
  phi0 = 0.2*twopi  ! 0.035*twopi has no resonance band with mL = 50
  kfloq = floquet_wavenumber(phi0)  ! phi0/2./2.**0.5 ! in units of m
  wn = floor(len_s*kfloq/twopi)
  print*,"w_n = ",wn
  
  !phi0 = 0.2*twopi, wn = 6; len_2=50, n=256 gives cool behaviour

  phi0 = 0.2*twopi; wn = 3
  
  call initialise_preheating_sine_wave(fld,phi0,1.e-4,wn)
  call time_evolve(twopi/128./(2.*nu**0.5), twopi*100./(2.*nu**0.5), twopi/16./(2.*nu**0.5))
#endif
  
#ifdef PREHEAT
  call initialise_model_symmetric(2, 1.e-1, 0._dl, 0._dl, 1000.)
  call setup_w_model(2, n, len_s/(2.*nu**0.5)/dble(n))
  fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  tcur => yvec(2*nLat*nFld+1)

  phi0 = 0.2*twopi
  kfloq = floquet_wavenumber(phi0)  ! phi0/2./2.**0.5
  print*,"floq index = ",kfloq*len_s/twopi

  nSamp = 1

  do i=1,nSamp
     call initialise_mean_preheating(fld, phi0)
     call initialise_fluctuations_white_long(fld, 1.e8, n/2) ! Check this (and convert to new version)
  
     call time_evolve(twopi/128./(2.*nu**0.5)/4., twopi*50./(2.*nu**0.5), twopi/16./(2.*nu**0.5))
  enddo
#endif
  
contains

  !>@brief Returns up limit of Floquet band in units of m, for sine-Gordon model
  real(dl) function floquet_wavenumber(phi0) result(kfloq)
    real(dl), intent(in) :: phi0

    kfloq = 0.5_dl*sqrt(0.5_dl)*phi0
  end function floquet_wavenumber
  
  subroutine initialise_false_vacuum(fld, spec_params, rho_bg)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: spec_params
    real(dl), intent(in) :: rho_bg

    fld = 0._dl
    call find_homogeneous_fv(fld, rho_bg)  ! Change the name of this
    yvec(4*nLat+1) = 0._dl
    call initialise_fluctuations(fld, spec_params)
  end subroutine initialise_false_vacuum

  ! Add model parameters in here
  !>@brief
  !> Find the false vacuum minima around which we will evolve our field fluctuations initially
  subroutine find_homogeneous_fv(fld,rho_bg)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: rho_bg
    real(dl) :: theta, sig

    sig = -1._dl
    ! This is just for the symmetric case
    fld(:,1,1) =  sqrt(rho_bg); fld(:,2,1) = 0._dl
    fld(:,1,2) = sig*sqrt(rho_bg); fld(:,2,2) = 0._dl

    ! Here is the result for nonsymmetric case
    ! Leading approximation
    theta = 0.125_dl*twopi + 0.25_dl*dg/(1._dl-gc - nu)
    theta = find_max_angle(nu,dg,gc)
    fld(:,1,1) =  (2._dl*rho)**0.5*cos(theta); fld(:,2,1) = 0._dl
    fld(:,1,2) = -(2._dl*rho)**0.5*sin(theta); fld(:,2,2) = 0._dl
  end subroutine find_homogeneous_fv

  
  ! Add input specifying the mean field
  subroutine sample_ics(nSamp, spec_params, nLat, nFld)
    integer, intent(in) :: nSamp
    type(SpecParams), intent(in) :: spec_params
    integer, intent(in) :: nLat, nFld
    
    integer :: o, i, j
    real(dl), dimension(1:nLat,1:2,1:nFld) :: fld_loc

    open(unit=newunit(o),file='initial_conditions.bin', access='stream')
    
    do i=1,nSamp
       fld_loc = 0._dl
       call initialise_fluctuations(fld_loc, spec_params)
       fld_loc(:,1,1) = fld_loc(:,1,1) + 1.  ! Change the mean to a parameter
       fld_loc(:,1,2) = fld_loc(:,1,2) - 1.  ! Change the mean to a passable parameter
       
       write(o) fld_loc
    enddo
    close(o)
  end subroutine sample_ics

  subroutine time_evolve_ensemble(nSamp, nLat, nFld, spec_params)
    integer, intent(in) :: nSamp
    integer, intent(in) :: nLat, nFld
    type(SpecParams) :: spec_params

    integer :: i

    ! Lots of nonlocality to kill in here
    do i=1,nSamp
       fld = 0._dl
       !call initialise_fluctuations_white_long(fld, rho, nCut=nCut)
       !call initialise_fluctuations_bogoliubov_long(fld, rho, nCut=nCut)
       call initialise_fluctuations_kg(fld, spec_params)
       fld(:,1,1) = 1._dl + fld(:,1,1)
       fld(:,1,2) = -1._dl + fld(:,1,1)

       !call time_evolve()
    enddo
  end subroutine time_evolve_ensemble
  
  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup
  !> before starting the program
  subroutine setup_w_model(nf,nl,dx)
    integer, intent(in) :: nf,nl
    real(dl), intent(in) :: dx
    call initialise_lattice(nf,nl,dx)
    call init_integrator(nVar)
  end subroutine setup_w_model

  ! Convert this to take a TimeStepper.  Or better, a Model object
  subroutine time_evolve(dt,tend,dtout)
    real(dl), intent(in) :: dt
    real(dl), intent(in) :: tend, dtout
    integer :: outsize, nums
    integer :: i,j

    if (dt > dx**2/twopi) print*,"Warning, potentially undersampling the Nyquist mode"

    outsize = floor(dtout/dt)
    nums = floor(tend/dt)/outsize
    
    dt_ = dt; dtout_ = dt*outsize ! Why the hell is this nonlocality here?  For output file, fix it!!!!
    call output_fields_binary(fld)
    call output_log_file(fld, 0.)
    
    tcur = 0._dl
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_log_file(fld, dt*i*outsize)
       call output_fields_binary(fld)
    enddo
  end subroutine time_evolve
  
  !>@brief
  !> Compute the desired time step dt adaptively using precomputed conditions.
  !>
  !>@todo
  !>@arg Write this.
  real(dl) function get_dt(dx,alph,w,wfrac) result(dt)
    real(dl), intent(in) :: dx, w
    real(dl), optional, intent(in) :: alph, wfrac
    real(dl) :: a_t, w_t
    real(dl) :: dt_cft, dt_w
    ! Need to satisfy:
    ! - Courant Condition
    ! - Resolve short time oscillator
    ! - Characteristic oscillation frequency

    a_t = 1._dl; w_t = 1._dl/16._dl
    if (present(alph)) a_t = alph
    if (present(wfrac)) w_t = wfrac
    
    dt_w = wfrac * twopi/w
    dt_cft = dx**2 / alph ! Fix this to take actual dispersion relationship
    dt = min(dt_w,dt_cft)
  end function get_dt

  !>@brief
  !> Initialise mean field with global relative phase phi0
  subroutine initialise_mean_preheating(fld, phi0, rho_norm)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: phi0
    real(dl), intent(in), optional :: rho_norm

    real(dl) :: rho_

    rho_ = 1._dl; if (present(rho_norm)) rho_ = rho_norm
    
    fld(:,1,1) = rho_; fld(:,2,1) = 0._dl
    fld(:,1,2) = rho_*cos(phi0); fld(:,2,2) = rho_*sin(phi0)
  end subroutine initialise_mean_preheating
  
  ! Subroutine to initialize fluctuations around constrained IR field
  ! This is copy/paste from above, actually need to rewrite it
  ! Move this into the fluctuations module
#ifdef CONSTRAINED
  subroutine initialise_constrained_fluctuations_IR(fld)
    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j

    spec = 0._dl
    do i=2,nLat/2
       spec(i) = (0.5_dl)**0.5 / sqrt(len*rho)
    enddo
    do i = 1,nFld; do j=1,2  ! i and j ordering is wrong
       call generate_1dGRF(df,spec(1:128),.false.)
       fld(:,i,j) = fld(:,i,j) + df
    enddo; enddo
  end subroutine initialise_constrained_fluctuations_IR
#endif

  subroutine initialise_preheating_sine_wave(fld,phi0,amp,wn)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: phi0
    real(dl), intent(in) :: amp
    integer, intent(in) :: wn

    real(dl), parameter :: rho_ave = 1._dl
    integer :: i

    do i=1,nLat
       fld(i,1,1) = sqrt(rho_ave); fld(i,2,1) = 0._dl
       fld(i,1,2) = sqrt(rho_ave)*cos(phi0+amp*sin(twopi*wn/len*(i-1)*dx)); fld(i,2,2) = sqrt(rho_ave)*sin(phi0+amp*sin(twopi*wn/len*(i-1)*dx))
    enddo
  end subroutine initialise_preheating_sine_wave

  subroutine initialise_preheating_density_wave(fld,phi0,amp,wn)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: phi0
    real(dl), intent(in) :: amp
    integer, intent(in) :: wn

    real(dl), dimension(1:nlat) :: drho_wave
    real(dl), parameter :: rho_ave = 1._dl
    integer :: i

    do i=1,nLat
       drho_wave = amp*sin(twopi*wn/len*(i-1)*dx)
    enddo

    fld(:,1,1) = sqrt(rho_ave*(1._dl-drho_wave)); fld(:,2,1) = 0._dl
    fld(:,1,2) = sqrt(rho_ave*(1._dl+drho_wave))*cos(phi0); fld(:,2,2) = sqrt(rho_ave*(1._dl+drho_wave))*sin(phi0) 
  end subroutine initialise_preheating_density_wave
  
  subroutine initialise_fields_sine(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), parameter :: rho_ave = 1._dl
    integer :: i; real(dl) :: dth, theta
    
    call initialise_mean_fields(fld,rho_ave)
    yvec(4*nLat+1) = 0._dl  ! turn this into tcur

    dth = twopi /dble(nLat)
    fld(:,1,1) = rho_ave; fld(:,2,1) = 0._dl  ! this should be a square root
    do i=1,nLat
       theta = 0.5_dl*twopi + (i-1)*dth
       fld(i,1,2) = rho_ave*cos(theta); fld(i,2,2) = rho*sin(theta)
    enddo
  end subroutine initialise_fields_sine

! Refactored into a separate file
#ifdef OUTPUT
!!!!!!!!!!!!!!!!!!!!
! Output Subroutines
!!!!!!!!!!!!!!!!!!!!

  subroutine output_log_file(fld,t,fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    real(dl), intent(in) :: t
    character(80), intent(in), optional :: fName

    character(80) :: fn
    integer, save :: oFile
    real(dl) :: en
    logical :: o
    
    fn = 'log.out'
    if (present(fName)) fn = trim(fName)
    
    inquire(file=trim(fn),opened=o)
    if (.not.o) then
       inquire(file=trim(fn),opened=o); if (o) close(oFile)
       open(unit=newunit(oFile), file=trim(fn))
       write(oFile,*) "# Lattice Parameters"
       write(oFile,*) "# n = ",nLat," dx = ",dx
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_
       write(oFile,*) "# Model Parameters"
       write(oFile,*) "# nu = ",nu," g = ",gs, " w = ",omega
       write(oFile,*) "# rho = ", rho
       write(oFile,*) "# delta = ", del
       write(oFile,*) "# t_{heal}   2\sqrt{nu}t_{heal} "
    endif
    write(oFile,*) t, 2.*sqrt(nu)*t
  end subroutine output_log_file
#endif
  
end program Gross_Pitaevskii_1d
