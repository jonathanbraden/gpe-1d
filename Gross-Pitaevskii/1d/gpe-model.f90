!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define FOURIER T

module gpe_model
  use constants, only : dl, twopi
#ifdef FOURIER
  use fftw3
#endif
  implicit none
  
  ! Properties of simulation grid
  real(dl) :: dx, dk, len
  integer :: nFld, nLat
  integer :: nVar
  real(dl), dimension(:), allocatable, target :: yvec 
#ifdef FOURIER
  type(transformPair1D) :: tPair 
#endif

  real(dl) :: omega
  ! Temporary parameter list.  Need to generalize.
  real(dl) :: gs, gc, dg
  real(dl) :: del, nu 
  real(dl) :: mu
  real(dl), allocatable :: lam(:,:)

  ! User specified model parameters.  Others are calculated
  type ModelParams
     integer :: nFld
     real(dl) :: nu, omega, del
     real(dl) :: gs, gc, dg
     real(dl) :: mu
     real(dl) :: rho_bg
  end type ModelParams

  ! Store parameters for the driving oscillator
  type NuParams
     real(dl) :: nu0, delta, omega
  end type NuParams
  
contains

  !>@brief
  !> Allocate storage for the lattice and current time
  !> Setup FFT for computing derivatives.
  subroutine create_lattice(nf,nl,dx_)
    integer, intent(in) :: nf, nl
    real(dl), intent(in) :: dx_

    nFld = nf; nLat = nl
    nVar = 2*nFld*nLat+1; allocate( yvec(1:nVar) )
    call initialize_transform_1d(tPair,nLat)
    
    dx = dx_; len=nLat*dx; dk = twopi/len
  end subroutine create_lattice

  function get_model_parameters() result(params)
    type(ModelParams) :: params

    params%nFld = nFld
    params%nu = nu; params%omega = omega; params%del = del
    params%gs = gs; params%gc = gc; params%dg = dg
    params%mu = mu
  end function get_model_parameters

  subroutine initialise_model_params(params)
    type(ModelParams), intent(in) :: params

    nFld = params%nFld
    nu = params%nu; omega = params%omega; del = params%del
    gs = params%gs; gc = params%gc; dg = params%dg
    mu = params%mu
  end subroutine initialise_model_params
    
  subroutine initialise_model_symmetric(nf,nu_,del_,om_)
    integer, intent(in) :: nf
    real(dl), intent(in) :: nu_, del_, om_

    nu = nu_
    gs=1._dl; gc = 0._dl; dg = 0._dl
    
    del = del_*(nu/2._dl)**0.5  
    omega = om_*2._dl*nu**0.5  ! conversion from scalar to healing length units
    !rho = rho_*2._dl*nu**0.5   ! conversion from scalar to healing length units

    ! Are my omega and rho conversions correct?
    call set_derived_params(nf)
  end subroutine initialise_model_symmetric
  
  subroutine initialise_model(nf,lv,nu_,gs_,gc_,dg_,om_)
    integer, intent(in) :: nf
    real(dl), intent(in) :: lv  ! Lambda parameter in scalar field potential
    real(dl), intent(in) :: nu_
    real(dl), intent(in) :: gs_, gc_, dg_
    real(dl), intent(in) :: om_

    gs = gs_; gc = gc_; dg = dg_
    nu = nu_

    del = (nu/2._dl)**0.5 * lv
    omega = om_*2._dl*nu**0.5
    call set_derived_params(nf)
  end subroutine initialise_model

  subroutine set_derived_params(nf)
    integer, intent(in) :: nf
    
    mu = 1._dl + nu  ! check if this is a plus or minus
    allocate(lam(1:nf,1:nf))
    lam = reshape( (/ gs+0.5_dl*dg,gc,gc,gs-0.5_dl*dg /), [nf,nf] )
  end subroutine set_derived_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conversions between dimensionless unit systems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief
  !> Convert from units of healing length and time (connected so program sound speed is sound speed of relative phase) to Drummonds units based on effective mass of scalar
  subroutine units_healing_to_scalar(nu_b)
    real(dl), intent(in) :: nu_b
    !omega_s = omega_h/(2._dl*nu_b**0.5)
    !len_s = len_h*(2._dl*nu_b)**0.5  ! check this one
    !rho_s = 0.5_dl*rho_h/(nu_b)**0.5 ! check this one
    !lam_s = del_h*(2._dl/nu_b)**0.5
  end subroutine units_healing_to_scalar
 
end module gpe_model
