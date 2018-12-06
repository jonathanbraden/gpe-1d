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
  
  ! GPE parameters
  integer :: nFld, nLat
  integer :: nVar
  real(dl), dimension(:), allocatable, target :: yvec  ! should remove target eventually

  ! Grid parameters
  real(dl) :: dx, dk, len
#ifdef FOURIER
  type(transformPair1D) :: tPair  ! needed in equation of motion, which now needs to ref this module
#endif

  real(dl) :: omega
  ! Temporary parameter list.  Need to generalize.
  real(dl) :: gs, gc, dg
  real(dl) :: del, nu, rho
  real(dl) :: mu
  real(dl), allocatable :: lam(:,:)

  ! User specified model parameters.  Others are calculated
  type mod_param
     integer :: nf
     real(dl) :: nu
     real(dl) :: omega
     real(dl) :: del, rho
     real(dl) :: gs, gc, dg
  end type mod_param
  
contains
  
  subroutine initialise_model_symmetric(nf,nu_,del_,om_,rho_)
    integer, intent(in) :: nf
    real(dl), intent(in) :: nu_, del_, rho_, om_

    nu = nu_
    gs=1._dl; gc = 0._dl; dg = 0._dl
    
    del = del_*(nu/2._dl)**0.5
    omega = om_*2._dl*nu**0.5  ! conversion from scalar to healing length units
    rho = rho_*2._dl*nu**0.5   ! conversion from scalar to healing length units

    call set_derived_params(nf)
  end subroutine initialise_model_symmetric
  
  subroutine initialise_model(nf,lv,nu_,gs_,gc_,dg_,om_,rho_)
    integer, intent(in) :: nf
    real(dl), intent(in) :: lv  ! Lambda parameter in scalar field potential
    real(dl), intent(in) :: nu_
    real(dl), intent(in) :: gs_, gc_, dg_
    real(dl), intent(in) :: om_,rho_

    gs = gs_; gc = gc_; dg = dg_
    nu = nu_

    del = (nu/2._dl)**0.5 * lv
    omega = om_*2._dl*nu**0.5
    rho = rho_*2._dl*nu**0.5
    call set_derived_params(nf)
  end subroutine initialise_model

  subroutine set_derived_params(nf)
    integer, intent(in) :: nf
    mu = 1._dl + nu  ! check if this is a plus or minus
    allocate(lam(1:nf,1:nf))
    lam = reshape( (/ gs+0.5_dl*dg,gc,gc,gs-0.5_dl*dg /), [nf,nf] )
  end subroutine set_derived_params
  
  subroutine initialise_lattice(nf,nl,dx_)
    integer, intent(in) :: nf, nl
    real(dl), intent(in) :: dx_

    nFld = nf; nLat = nl
    nVar = 2*nFld*nLat+1; allocate( yvec(1:nVar) )
    print*,nVar,nFld
    call initialize_transform_1d(tPair,nLat)
    
    dx = dx_; len=nLat*dx; dk = twopi/len
  end subroutine initialise_lattice

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

  !>@brief
  !> Convert from Drummond's units to healing length and time
  subroutine units_scalar_to_healing(nu_b)
    real(dl), intent(in) :: nu_b
    !omega_h = omega_2 * 2._dl*nu_b**0.5
  end subroutine units_scalar_to_healing
  
end module gpe_model
