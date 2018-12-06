program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use gaussianRandomField
  use eom
  use integrator
  implicit none
  real(dl), dimension(:,:,:), pointer :: fld
  real(dl), pointer :: tcur
  real(dl) :: dtout_, dt_

  integer :: n
  ! Time stepping parameters
  integer :: ds, w_samp

  ! Fixing parameters for now
  real(dl), parameter :: omega_s=50._dl, rho_s=1000._dl, del_s=1._dl+0.4_dl
  real(dl), parameter :: len_s = 50._dl
  
  n=256
  ds = 1; w_samp=16
  
  call initialise_model(2,del_s,2.e-3,1._dl,0._dl,0._dl,omega_s,rho_s)
  call setup_w_model(2,n,len_s/(2.*(nu)**0.5)/dble(n))
  fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  tcur => yvec(2*nLat*nFld+1)
  call initialise_fields_rand(fld)
  
  call time_evolve((1._dl/dble(w_samp))*twopi/omega,100._dl,dble(ds)/dble(w_samp)*twopi/omega)
  
contains

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

    if (dt > dx**2) print*,"Warning, violating Courant condition"

    outsize = floor(dtout/dt)
    nums = floor(tend/dt)/outsize
    
    dt_ = dt; dtout_ = dt*outsize ! Why the hell is this nonlocality here?  For output file, fix it!!!!
    call output_fields(fld,0.)

    tcur = 0._dl
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_fields(fld,tcur)
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
  !> Find the false vacuum minima around which we will evolve our field fluctuations initially
  subroutine initialise_mean_fields(fld,rho)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: rho
    real(dl) :: theta, sig

    sig = -1._dl
    ! This is just for the symmetric case
    fld(:,1,1) =  rho; fld(:,2,1) = 0._dl
    fld(:,1,2) = -rho; fld(:,2,2) = 0._dl
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
  
  subroutine initialise_fluctuations(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld

    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j
    
    spec = 0._dl
    do i=2,nLat/2
       spec(i) = (0.5_dl)**0.5 / (sqrt(len*rho))  ! This is the only place I need rho
    enddo
    do i = 1,nFld; do j=1,2
       call generate_1dGRF(df,spec(1:128),.false.)
       fld(:,i,j) = fld(:,i,j) + df
    enddo; enddo
  end subroutine initialise_fluctuations

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
    do i = 1,nFld; do j=1,2
       call generate_1dGRF(df,spec(1:128),.false.)
       fld(:,i,j) = fld(:,i,j) + df
    enddo; enddo
  end subroutine initialise_constrained_fluctuations_IR
#endif

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
    real(dl), parameter :: rho = 1._dl
    integer :: i; real(dl) :: dth, theta
    
    call initialise_mean_fields(fld,rho)
    yvec(4*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialise_fluctuations(fld)
  end subroutine initialise_fields_rand
    
  subroutine output_fields(fld,t,fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    real(dl), intent(in) :: t
    character(80), intent(in), optional :: fName
    logical :: o; integer :: i
    character(80) :: fn

    integer, parameter :: oFile = 99  ! remove this ugliness
    
    fn = 'fields.dat'
    if (present(fName)) fn = trim(fName)
    
    inquire(file=trim(fn),opened=o)
    if (.not.o) then
       inquire(file=trim(fn),opened=o); if (o) close(oFile)
       open(unit=oFile,file=trim(fn))
       write(oFile,*) "# Lattice Parameters"
       write(oFile,*) "# n = ",nLat," dx = ",dx
       write(oFile,*) "# Time Stepping parameters"
       write(oFile,*) "# dt = ",dt_, " dt_out = ",dtout_
       write(oFile,*) "# Model Parameters"
       write(oFile,*) "# nu = ",nu," g = ",gs, " w = ",omega
       write(oFile,*) "# rho = ", rho
       write(oFile,*) "# delta = ", del
    endif
       
    do i=1,nLat
       write(oFile,*) fld(i,:,:)
    enddo
    write(oFile,*)
    
  end subroutine output_fields
  
end program Gross_Pitaevskii_1d
