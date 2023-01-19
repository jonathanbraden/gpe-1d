program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use utils, only : newunit
  use gaussianRandomField
  use eom
  use output
  use integrator
  implicit none
  real(dl), dimension(:,:,:), pointer :: fld
  real(dl), pointer :: tcur
  real(dl) :: dtout_, dt_

  integer :: n
  integer :: ds, w_samp

  ! Fixing parameters for now
  real(dl), parameter :: omega_s = 50._dl, rho_s = 1000._dl, del_s = 1._dl+0.4_dl
  real(dl), parameter :: len_s = 100._dl

  real(dl) :: phi0, kfloq
  integer :: wn

  integer :: i, nSamp

  real(dl) :: om_alex, rho_alex, nu_alex, del_alex, phi_init, len_alex
  
  n=512
  ds = 1; w_samp=16

  ! Old stuff for false vacuum
  !call initialise_model(2, del_s, 2.e-3, 1._dl, 0.5_dl, 0.5_dl, omega_s, rho_s)
  !call setup_w_model( 2, n, len_s/(2.*(nu)**0.5)/dble(n) )
  !fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  !tcur => yvec(2*nLat*nFld+1)
  !call initialise_fields_rand(fld)
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

  
  !call sample_ics(1000, n, 2, n/2, rho_alex)
  
  ! Initialise the fluctuations
  fld = 0._dl
  !call initialise_fluctuations_white(fld, rho_alex, nCut=n/2)
  !call initialise_fluctuations_bogoliubov(fld, rho_alex, nCut=n/2)
  call initialise_fluctuations_kg(fld, rho, nCut=n/2, modes=(/.true.,.false./))

  ! Add displaced mean field
  fld(:,1,1) = fld(:,1,1) + 1.
  fld(:,1,2) = fld(:,1,2) - 1.
  !fld(:,1,1) = cos(0.5*phi_init) + fld(:,1,1); fld(:,2,1) = -sin(0.5*phi_init) + fld(:,2,1)
  !fld(:,1,2) = -cos(0.5*phi_init) + fld(:,1,2); fld(:,2,2) = -sin(0.5*phi_init) + fld(:,2,2)

  !print*,"lam is ",lam," del = ",del," nu = ",nu, "mu = ",mu
  
  call time_evolve( (twopi/om_alex)/dble(w_samp), 2.*45.465, (twopi/om_alex)/dble(2.) )
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
  kfloq = phi0/2./2.**0.5 ! in units of m
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
  kfloq = phi0/2./2.**0.5
  print*,"floq index = ",kfloq*len_s/twopi

  nSamp = 1

  do i=1,nSamp
     call initialise_mean_preheating(fld, phi0)
     call initialise_fluctuations_white(fld, 1.e8)
  
     call time_evolve(twopi/128./(2.*nu**0.5)/4., twopi*50./(2.*nu**0.5), twopi/16./(2.*nu**0.5))
  enddo
#endif
  
contains

  subroutine sample_ics(nSamp, nLat, nFld, nCut, rho)
    integer, intent(in) :: nSamp
    integer, intent(in) :: nLat, nFld, nCut
    real(dl), intent(in) :: rho
    
    integer :: o, i, j
    real(dl), dimension(1:nLat,1:2,1:nFld) :: fld_loc

    ! Add input specifying the mean field
    o=70
    open(unit=o,file='initial_conditions.bin', access='stream')
    
    do i=1,nSamp
       fld_loc = 0._dl
       !call initialise_fluctuations_white(fld_loc, rho, nCut=nCut)
       !call initialise_fluctuations_bogoliubov(fld_loc, rho, nCut=nCut)
       call initialise_fluctuations_kg( fld_loc, rho, nCut=nCut, modes=(/.true.,.false./) )
       fld_loc(:,1,1) = fld_loc(:,1,1) + 1.
       fld_loc(:,1,2) = fld_loc(:,1,2) - 1.
       write(o) fld_loc
    enddo
    close(o)
  end subroutine sample_ics

  subroutine time_evolve_ensemble(nSamp, nLat, nFld, nCut, rho)
    integer, intent(in) :: nSamp
    integer, intent(in) :: nLat, nFld, nCut
    real(dl), intent(in) :: rho

    integer :: i

    
    do i=1,nSamp
       fld = 0._dl
       !call initialise_fluctuations(fld, rho, nCut=nCut)
       !call initialise_fluctuations_bogoliubov(fld, rho, nCut=nCut)
       call initialise_fluctuations_kg(fld, rho, nCut=nCut)
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

  ! Moved to a module
  real(dl) function dispersion_total(k) result(om)
    real(dl), intent(in) :: k

    om = k*sqrt(1. + k**2/4._dl)
  end function dispersion_total

  ! Moved to a module
  real(dl) function dispersion_rel(k,nu) result(om)
    real(dl), intent(in) :: k, nu
    om = 2._dl*sqrt( (nu+0.25*k**2) * (1._dl+nu+0.25*k**2) )
  end function dispersion_rel
  
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
  
  !>@brief
  !> Find the false vacuum minima around which we will evolve our field fluctuations initially
  subroutine initialise_mean_fields(fld,rho)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: rho
    real(dl) :: theta, sig

    sig = -1._dl
    ! This is just for the symmetric case
    fld(:,1,1) =  rho**0.5; fld(:,2,1) = 0._dl
    fld(:,1,2) = -rho**0.5; fld(:,2,2) = 0._dl

    ! Here is the result for nonsymmetric case
    ! Leading approximation
    theta = 0.125_dl*twopi + 0.25_dl*dg/(1._dl-gc - nu)
    theta = find_max_angle(nu,dg,gc)
    fld(:,1,1) =  (2._dl*rho)**0.5*cos(theta); fld(:,2,1) = 0._dl
    fld(:,1,2) = -(2._dl*rho)**0.5*sin(theta); fld(:,2,2) = 0._dl
    print*,"Homogeneous ICs are ",fld(1,1,1),fld(1,1,2)
  end subroutine initialise_mean_fields

  real(dl) function find_max_angle(nu,dg,gc) result(theta)
    real(dl), intent(in) :: nu, dg, gc
    real(dl) :: nu_bg, dg_bg, dth; integer :: l
    integer, parameter :: maxit = 16; real(dl), parameter :: eps = 1.e-16
    print*,"in the function"
    theta = 0.125_dl*twopi + 0.25_dl*dg/(1._dl-gc - nu)
    nu_bg = nu / (1._dl-gc); dg_bg = dg / (1._dl-gc)
    print*,"entering loop"
    do l = 1,maxit
       dth = - dvTheta(theta,nu_bg,dg_bg) / d2vTheta(theta,nu_bg,dg_bg)
       theta = theta + dth
       if (abs(dth) < eps) exit
    end do
    print*,"Converged in ",l," iterations to ",theta
    if (l==maxit) print*,"Newton method failed to find a background"
  end function find_max_angle

  real(dl) function theta_guess(nu,dg,gc) result(theta)
    real(dl), intent(in) :: nu,dg,gc
    theta = 0.125_dl*twopi + 0.25_dl*dg/(1._dl-gc-nu)
  end function theta_guess
    
  real(dl) function dvTheta(theta, nu_bg, dg_bg) result(dv)
    real(dl), intent(in) :: theta, nu_bg, dg_bg
    real(dl), parameter :: sig = -1._dl
    dv = cos(2._dl*theta)*sin(2._dl*theta) + sig*nu_bg*cos(2.*theta) + 0.5*dg_bg*sin(2.*theta)
  end function dvTheta

  real(dl) function d2vTheta(theta, nu_bg, dg_bg) result(d2v)
    real(dl), intent(in) :: theta, nu_bg, dg_bg
    real(dl), parameter :: sig = -1._dl
    d2v = 2._dl*(cos(2._dl*theta)**2 - sin(2._dl*theta)**2) - 2.*sig*nu_bg*sin(2.*theta) + dg_bg*cos(2.*theta)
  end function d2vTheta
  
  subroutine initialise_fluctuations_white(fld, rho, nCut_)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: rho
    integer, intent(in), optional :: nCut_

    real(dl) :: df(1:size(fld,dim=1)), spec(1:size(fld,dim=1)/2+1)
    integer :: i,j
    integer :: n, num_fld
    integer :: nCut

    n = size(fld,dim=1)
    num_fld = size(fld,dim=3)
    nCut = n/2 ;  if (present(nCut_)) nCut = nCut_
    
    spec = 0._dl
    do i=2,n/2
       spec(i) = 1._dl / (sqrt(2._dl*len*rho))  ! This is the only place I need rho
    enddo
    do i = 1,num_fld; do j=1,2
       call generate_1dGRF(df,spec(1:nCut),.false.)
       fld(:,j,i) = fld(:,j,i) + df
    enddo; enddo
  end subroutine initialise_fluctuations_white

  subroutine initialise_fluctuations_kg(fld, rho, nCut, modes)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: rho
    integer, intent(in) :: nCut
    logical, intent(in), dimension(2), optional :: modes

    real(dl), dimension(1:size(fld,dim=1),1:2) :: df_rel, df_tot
    real(dl) :: spec(1:size(fld,dim=1)/2+1)
    real(dl) :: norm, dk, keff
    real(dl) :: m2eff, lameff
    integer :: i,j
    logical, dimension(2) :: modes_
    
    modes_ = (/ .true., .true. /);  if (present(modes)) modes_ = modes
    
    lameff = del*sqrt(2._dl/nu)
    m2eff = 4.*nu*(lameff**2-1._dl)
    norm = 1._dl / sqrt(2._dl*len*rho)
    dk = twopi/len

    df_rel = 0._dl
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl/sqrt(keff)
    enddo
    spec = spec * norm
    do i = 1,nFld; do j=1,2
       if (modes(2)) call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo; enddo

    df_tot = 0.
    spec = 0._dl
    do i=1,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl / (keff**2+m2eff)**0.25 
    enddo
    spec = spec * norm
    do i=1,nFld; do j=1,2
       if (modes(1)) call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo; enddo

    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo
  end subroutine initialise_fluctuations_kg
    
  subroutine initialise_fluctuations_bogoliubov(fld, rho, nCut)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: rho
    integer, intent(in), optional :: nCut

    real(dl) :: df_rel(1:nLat,1:2), df_tot(1:nLat,1:2), spec(1:nLat/2+1)
    integer :: i,j
    real(dl) :: norm, dk, keff
    real(dl) :: lameff, nu_

    nu_ = nu  ! Fix this nonlocality
    lameff = del*sqrt(2._dl/nu_)
    norm = 1._dl / sqrt(2._dl*len*rho)
    dk = twopi / len   ! Fix this nonlocality
    
    ! Relative modes
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = sqrt( (keff**2+2._dl)/(keff*sqrt(keff**2+4._dl)) )
    enddo
    spec = spec * norm
    do j=1,2
       call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo
    
    ! Now get total modes
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = sqrt( (keff**2 + 2._dl - 4._dl*nu_)  / sqrt(keff**2+4._dl*nu_*(lameff**2-1._dl)) / sqrt(keff**2+4._dl-4._dl*nu_*(lameff**2+1._dl)) )
    enddo
    spec = spec * norm
    do j=1,2
       call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo
    
    ! Now remix the fields and add them to the background
    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo
    
  end subroutine initialise_fluctuations_bogoliubov
  
  ! Subroutine to initialize fluctuations around constrained IR field
  ! This is copy/paste from above, actually need to rewrite it
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
  
  subroutine initialise_fields_rand(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), parameter :: rho_bg = 1._dl
    integer :: i; real(dl) :: dth, theta
    real(dl) :: rho_fluc

    rho_fluc = rho ! Fix this horrible nonlocality
    call initialise_mean_fields(fld,rho_bg)
    yvec(4*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialise_fluctuations_white(fld, rho_fluc)  ! Fix this horrible nonlocality
  end subroutine initialise_fields_rand

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
  
  subroutine output_fields_binary(fld,fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    character(80), intent(in), optional :: fName

    logical :: o
    character(80) :: fn
    integer, save :: oFile
    
    fn = 'fields.bin'
    if (present(fName)) fn = trim(fName)

    inquire(file=trim(fn), opened=o)
    if (.not.o) open(unit=newunit(oFile), file=fn, access='stream')

    write(oFile) fld
  end subroutine output_fields_binary
#endif
  
end program Gross_Pitaevskii_1d
