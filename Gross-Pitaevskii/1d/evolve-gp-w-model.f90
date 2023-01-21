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

  real(dl), dimension(:,:,:), pointer :: fld
  real(dl), pointer :: tcur
  type(SpecParams) :: fluc_params
  type(TimeStepper) :: time_stepper

  integer :: i ! Automatic creation of array
  
  ! These are input parameters used for readability
  integer :: n

  ! Fixing parameters for now
  real(dl), parameter :: omega_s = 50._dl, del_s = 1._dl + 0.4_dl !rho_s = 1000._dl
  real(dl), parameter :: len_s = 100._dl

  ! Used for preheating code
  real(dl) :: phi0, kfloq
  integer :: wn
  integer :: nSamp

  ! Move params somewhere else
  real(dl) :: om_gpe, nu_gpe, del_gpe, len_gpe
  real(dl) :: rho_alex, phi_init

  integer :: ds, w_samp
  
  
  ! Old stuff for false vacuum
  !n=512
  !ds = 1; w_samp=16
  !call initialise_model(2, del_s, 2.e-3, 1._dl, 0.5_dl, 0.5_dl, omega_s, rho_s)
  !call setup_simulation( 2, n, len_s/(2.*(nu)**0.5)/dble(n), fld, tcu )
  !call initialise_false_vacuum(fld)
  !call time_evolve((1._dl/dble(w_samp))*twopi/omega,100._dl,dble(ds)/dble(w_samp)*twopi/omega)

  n = 512
  nu_gpe = 0.01
  del_gpe = 1.2_dl
  om_gpe = 64./2./sqrt(nu_gpe)
  
  len_gpe = 50./2./sqrt(nu_gpe) !454.65
  
  w_samp = 16
  rho_alex = 43.99 !43.99 !43.99 ! This is only really needed for the ICs, not the model setup 
  phi_init = 0._dl  ! 0.01*twopi

  call setup_simulation( 2, n, len_gpe/n, fld, tcur )
  call initialise_model_symmetric( 2, nu_gpe, del_gpe, om_gpe ) ! fix rho, num-field nonlocality

  !fluc_params = make_spec_params( rho_alex, len_gpe, nu_gpe, del_gpe, 'BOGO', (/( .true., i=1,size(fld,dim=3))/), n/2, fv_=.true. )
  fluc_params = make_spec_params( rho_alex, len_gpe, nu_gpe, del_gpe, 'BOGO', (/ .true., .false. /), n/2, fv_=.true. )
  call print_spec_params(fluc_params)

  !!! Sample initial condtions
  nSamp = 1000
  call sample_ics(nSamp, fluc_params, n, 2)

  nSamp = 10
  do i=1,nSamp
     fld = 0._dl
     call initialise_fluctuations(fld, fluc_params)
     fld(:,1,1) = fld(:,1,1) + 1.; fld(:,1,2) = fld(:,1,2) - 1.
     !fld(:,1,1) = fld(:,1,1) + 1.; fld(:,1,2) = fld(:,1,2) + 1.
  !fld(:,1,1) = cos(0.5*phi_init) + fld(:,1,1); fld(:,2,1) = -sin(0.5*phi_init) + fld(:,2,1)
  !fld(:,1,2) = -cos(0.5*phi_init) + fld(:,1,2); fld(:,2,2) = -sin(0.5*phi_init) + fld(:,2,2)

     call set_time_steps_oscillator(time_stepper, om_gpe, w_samp, out_size=8*w_samp, t_final=len_gpe)  ! 45.456
     call print_time_stepper(time_stepper)
     call time_evolve_stepper(fld, time_stepper, verbose_=.false.)
  enddo
     
! With improved interface, rewrite this
#ifdef PREHEAT_SINGLE
  call initialise_model_symmetric( 2, 1.e-2, 0._dl, 0._dl )
  call setup_simulation(2, n, len_s/(2.*(nu)**0.5)/dble(n), fld, tcur)
  
  phi0 = 0.2*twopi  ! 0.035*twopi has no resonance band with mL = 50
  kfloq = floquet_wavenumber(phi0)  ! phi0/2./2.**0.5 ! in units of m
  wn = floor(len_s*kfloq/twopi)
  print*,"w_n = ",wn
  
  !phi0 = 0.2*twopi, wn = 6; len_2=50, n=256 gives cool behaviour
  phi0 = 0.2*twopi; wn = 3
  
  call initialise_preheating_sine_wave(fld,phi0,1.e-4,wn)
  call time_evolve(twopi/128./(2.*nu**0.5), twopi*100./(2.*nu**0.5), twopi/16./(2.*nu**0.5))
#endif

! With improved interface, rewrite this
#ifdef PREHEAT
  call initialise_model_symmetric(2, 1.e-1, 0._dl, 0._dl )  ! 1000 is the value for rho
  call setup_simulation(2, n, len_s/(2.*nu**0.5)/dble(n), fld, tcur)

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

  subroutine evolve_ensemble(nSamp, stepper, spec_params)
    integer, intent(in) :: nSamp
    type(TimeStepper), intent(inout) :: stepper
    type(SpecParams), intent(in) :: spec_params

    real(dl), dimension(1:size(fld,dim=3)) :: mean_fld
    integer :: i,j

    print*,"Running an ensemble of ",nSamp," realisations"
    call print_time_stepper(stepper)
    call print_spec_params(spec_params)

    mean_fld = 0._dl
    mean_fld(1) = 1._dl
    mean_fld(2) = -1._dl
    
    do i=1,nSamp
       fld = 0._dl
       call initialise_fluctuations(fld, spec_params)
       do j=1,size(fld,dim=3)
          fld(:,1,j) = fld(:,1,j) + mean_fld(j)
       enddo
       call time_evolve_stepper(fld, stepper, verbose_=.false.)
    enddo
  end subroutine evolve_ensemble
  
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
    fld(:,1,1) =  (2._dl*rho_bg)**0.5*cos(theta); fld(:,2,1) = 0._dl
    fld(:,1,2) = -(2._dl*rho_bg)**0.5*sin(theta); fld(:,2,2) = 0._dl
  end subroutine find_homogeneous_fv

  
  ! Add input specifying the mean field
  subroutine sample_ics(nSamp, spec_params, nLat, nFld)
    integer, intent(in) :: nSamp
    type(SpecParams), intent(in) :: spec_params
    integer, intent(in) :: nLat, nFld
    
    integer :: o, i
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
  
  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup
  !> before starting the program
  subroutine setup_simulation(nf,nl,dx, fld, tcur)
    integer, intent(in) :: nf,nl
    real(dl), intent(in) :: dx
    real(dl), dimension(:,:,:), pointer :: fld
    real(dl), pointer :: tcur
    
    call create_lattice(nf,nl,dx)
    call init_integrator(nVar)

    fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
    tcur => yvec(2*nLat*nFld+1)
  end subroutine setup_simulation
  
  ! Convert this to take a TimeStepper.  Or better, a Model object
  subroutine time_evolve(dt,tend,dtout)
    real(dl), intent(in) :: dt
    real(dl), intent(in) :: tend, dtout
    integer :: outsize, nums
    integer :: i,j

    if (dt > dx**2/twopi) print*,"Warning, potentially undersampling the Nyquist mode"

    outsize = floor(dtout/dt)
    nums = floor(tend/dt)/outsize
    
    call output_fields_binary(fld)
    call output_log_file(fld, 0., dt)
    
    tcur = 0._dl
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_log_file(fld, dt*i*outsize, dt)
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

  ! Do I use this or need it?  It's clearly got buggy behaviour with the rho
  subroutine initialise_fields_sine(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), parameter :: rho_bg = 1._dl
    integer :: i; real(dl) :: dth, theta
    
    call initialise_mean_fields(fld,rho_bg)
    yvec(4*nLat+1) = 0._dl  ! turn this into tcur

    dth = twopi /dble(nLat)
    fld(:,1,1) = rho_bg; fld(:,2,1) = 0._dl  ! this should be a square root
    do i=1,nLat
       theta = 0.5_dl*twopi + (i-1)*dth
       fld(i,1,2) = rho_bg*cos(theta); fld(i,2,2) = rho_bg*sin(theta)
    enddo
  end subroutine initialise_fields_sine
  
end program Gross_Pitaevskii_1d
