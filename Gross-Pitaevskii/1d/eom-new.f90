!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the variables and equations describing the model we are solving
!>
!> This module provides storage for the complex fields describing a collection
!> of Bose condensates, along with appropriate time evolution routines.
!> We assume the condensates obey the Gross-Pitaevskii equation
!>  \f[
!>    i\hbar \frac{\partial\psi_i}{\partial t} = -\frac{\hbar^2}{2m}\nabla^2\psi_i - \mu_i\psi_i + \sum_j g_{ij}\left|\psi_j\right|^2\psi_i - \sum_j\nu_{ij}psi_j
!>  \f]
!> with the coefficients \f$g_{ij},\mu_i\f$ and \f$\nu_{ij}\f$ viewed as model parameters
!> to be adjusted in the experimental setup.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define ADJUST T
#define FOURIER T
#define DISCRETE T
!#define FOURIER_DIFF T
#define DIFF T
#define FIELD_TYPE Field_Model
module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  implicit none

  integer, parameter :: nFld = 2, nLat = 512
  integer, parameter :: nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec

  ! Parameters for 2-field model in "scalar field" normalisation
  real(dl) :: nu_s
  real(dl), parameter :: omega_s = 50._dl*1., rho_s=1000._dl, del_s= 1._dl + 0.3_dl
  real(dl), parameter :: len_s=50._dl

  real(dl), parameter :: nu = 2.e-3, mu=1._dl+nu, rho = rho_s*2._dl*(nu)**0.5
  real(dl), parameter :: omega = omega_s*2._dl*nu**0.5, gc=0.1_dl, gs=1._dl, dg=0.43_dl
  real(dl), parameter :: del = (nu/2._dl)**0.5 * del_s
  real(dl), parameter :: delm = 0._dl  ! Mass difference
  
  real(dl), parameter, dimension(1:nFld) :: m = (/ 1._dl, 1._dl /)
  real(dl), dimension(1:nFld,1:nFld), parameter :: lam = reshape( (/ gs+0.5_dl*dg, gc, gc, gs-0.5_dl*dg /), [nFld,nFld] )
  
  real(dl), parameter :: len = len_s / (2.*(nu)**0.5), dx = len/dble(nLat), dk = twopi/len
  ! Drummond parameters
  real(dl) :: nu_scl, rho_scl, delm_scl, lam_scl, omega_scl 
  
#ifdef FOURIER
  type(transformPair1D) :: tPair
#endif
  
contains

#ifdef MODDED
  !>@brief
  !> Allocate storage space for the field lattice variables and parameters of the lattice
  subroutine create_field_lattice(this,n,nf,len)
    type(FIELD_TYPE), intent(out) :: this
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: len

    this%nx = n; this%nfld = nf
    this%len = len
    this%dx = len / dble(n); this%dk = twopi / len
    this%nvar = 2*nfld*n
    allocate(this%yvec(1:nvar))
  end subroutine create_field_lattice
#endif
  ! Rather than hardcoding parameters above, use this subroutine
  subroutine initialise_field(n)
    integer, intent(in) :: n
    !call initialise_transform_1d(tPair,n)
  end subroutine initialise_field
  
  !>@brief
  !> Converts the dimensionless coupling constants from those associated with the scalar potential
  !> \f[
  !>    \hbar\omega_p = 2\sqrt{\nu gn} \qquad c_p^2 = \frac{gn}{m} \qquad \kappa_p = c_p^{-1}\omega_p
  !> \f]
  !> to the "condensate" time units used in this code
  !> \f[ 
  !>   \hbar\omega_p = gn \qquad c_p^2 = \frac{gn}{m} \qquad \kappa_p = c_p^{-1}\omega_p
  !> \f]
  !>
  !> To Do: Rename this to convert from scalar potential
  !>  Will need, V_0 = w_0^2, lambda, and tilde(nu), etc
  subroutine convert_units_scalar_to_condensate(nu_b,l_d,w_d,rho_d,lam_d)
    real(dl), intent(in) :: nu_b, l_d, w_d, rho_d, lam_d
    real(dl) :: w_c, l_c, rho_c, del_c
    
    w_c = w_d * 2._dl*(nu_b)**0.5
    l_c = l_d / (2._dl*nu_b**0.5)
    rho_c = rho_d * 2._dl*(nu_b)**0.5
    del_c = (nu_b/2._dl)**0.5*lam_d
  end subroutine convert_units_scalar_to_condensate

  !>@brief
  !> Convert dimensionless units from those associated with the condensate to those associated with the scalar field.
  !> Performs the inverse operation to convert_units_scalar_to_condensate.
  !>
  !>@todo
  !> Write this
  subroutine convert_units_condensate_to_scalar(nu_b,del_c,rho_c,l_c,w_c)
    real(dl), intent(in) :: nu_b, del_c, rho_c, l_c, w_c
    real(dl) :: w_s, l_s, rho_s, lam_s

    w_s = 0.5_dl*w_c / nu_b**0.5
    l_s = l_c*(2._dl*nu_b)**0.5
    rho_s = 0.5_dl*rho_c/(nu_b)**0.5
    lam_s = del_c*(2._dl/nu_b)**0.5
  end subroutine convert_units_condensate_to_scalar
  
  subroutine set_model_params_2fld(g_i,dg_i,gc_i,nu_i)
    real(dl), intent(in), optional :: g_i, dg_i, gc_i, nu_i
    real(dl) :: g_t, dg_t, gc_t, nu_t

    g_t = 0._dl; dg_t = 0._dl; gc_t = 0._dl; nu_t = 0._dl
    ! Add all the appropriate if statements
  end subroutine set_model_params_2fld

  !>@brief
  !> Solve the background condensate densities in the 2-field model
  function background_densities(dg) result(n)
    real(dl), intent(in) :: dg
    real(dl), dimension(1:nFld) :: n
    n = 1._dl  ! Fix this to actually solve
  end function background_densities

! Derivatives of potential with respect to theta for computing mean densities by Newton's method
  function dVdTheta(x) result(dv)
    real(dl), intent(in) :: x
    real(dl) :: dv

    dv = cos(2._dl*x)*(sin(2._dl*x) + 4._dl*nu) + sin(2._dl*x)*dg
  end function dVdTheta

  function d2Vd2Theta(x) result(d2v)
    real(dl), intent(in) :: x
    real(dl) :: d2v
    
  end function d2Vd2Theta
  
#ifdef ADJUST
  !>@brief
  !> Set the parameters of the Bose-Condensates assuming the simplified situation of equal couplings
  subroutine set_model_params_simple(mu_i,g_i,gc_i,nu_i)
    real(dl), intent(in), optional :: mu_i,g_i,gc_i,nu_i
    integer :: i
    
    if ( present(mu_i) )       mu = mu_i
    if ( present(g_i) )  then; do i=1,nFld; lam(i,i) = g_i; enddo; endif
    if ( present(gc_i) ) then; endif  ! Finish this one
    if ( present(nu_i) ) then; nu = nu_i; do i=1,nFld; nu(i,i) = 0._dl; enddo; endif
  end subroutine set_model_params_simple
  
  subroutine set_model_params(mu_i,g_i,nu_i)
    real(dl), optional, intent(in) :: mu_i, g_i(1:nFld,1:nFld), nu_i(1:nFld,1:nFld)

    if (present(mu_i)) mu = mu_i
    if (present(g_i)) lam = g_i
    if (present(nu_i)) nu = nu_i
  end subroutine set_model_params
#endif

  subroutine set_drummond_params(nu_i,rho_i,lam_i, omega_i)
    real(dl), intent(in) :: nu_i, rho_i, lam_i, omega_i 
    nu_scl = nu_i; rho_scl = rho_i; lam_scl = lam_i; omega_scl = omega_i
  end subroutine set_drummond_params
  
#define R1 1:nLat
#define I1 (nLat+1):(2*nLat)
#define R2 (2*nLat+1):(3*nLat)
#define I2 (3*nLat+1):(4*nLat)

#define TIND 4*nLat+1
  !>@brief
  !> Compute the derivatives of the the vector of fields
  !>
  !> Compute the time derivatives of the real and imaginary parts of the
  !> Bose-Condensate field \f$\psi_i\f$
  !>  \f[
  !>    \psi_i^R =
  !>  \f]
  !>  \f[
  !>    \psi_i^I = 
  !>  \f]
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    real(dl) :: nueff
    real(dl), dimension(1:nLat,1:nFld) :: ysq
#ifdef DISCRETE
    real(dl), parameter :: lNorm = 0.5_dl/dx**2
#endif
    
    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2

    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    nueff = nu + del*omega*cos(omega*yc(TIND))
    
    ! These are ugly, add some default vectors so vectorisation can be done more easily
    yp(R1) = -mu*yc(I1) + ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(I1) ! -laplacian(I1)
    yp(I1) = mu*yc(R1)  - ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(R1) ! +laplacian(R1)
    yp(R2) = -mu*yc(I2) + ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(I2)  ! -laplacian(I2)
    yp(I2) = mu*yc(R2)  - ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(R2) ! + laplacian(R2)

    yp(R1) = yp(R1) - nueff*yc(I2)
    yp(I1) = yp(I1) + nueff*yc(R2)
    yp(R2) = yp(R2) - nueff*yc(I1)
    yp(I2) = yp(I2) + nueff*yc(R1)

    ! Fourier based differentiation
    ! Need to add Chebyshev in here
#ifdef DIFF
#ifdef DISCRETE
    ! Compute the endpoints
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    yp(1) = yp(1) - lNorm * ( yc(2*nlat) - 2._dl*yc(nlat+1) + yc(nlat+2) )
    yp(nlat) = yp(nlat) - lNorm * ( yc(2*nlat-1) - 2._dl*yc(2*nlat) + yc(nlat+1) )
    yp(3*nlat+1) = yp(3*nlat+1) + lNorm * ( yc(3*nlat) - 2._dl*yc(2*nlat+1) + yc(2*nlat+2) )
    yp(4*nlat) = yp(4*nlat) + lNorm * ( yc(3*nlat-1) - 2._dl*yc(3*nlat) + yc(2*nlat+1) )
    yp(2*nlat+1) = yp(2*nlat+1) - lNorm * ( yc(4*nlat) - 2._dl*yc(3*nlat+1) + yc(3*nlat+2) )
    yp(3*nlat) = yp(3*nlat) - lNorm * ( yc(4*nlat-1) - 2._dl*yc(4*nlat) + yc(3*nlat+1) )
    
    yp(nlat+2:2*nlat-1) =   yp(nlat+2:2*nlat-1)     + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
    yp(2:nlat-1) =          yp(2:nlat-1)   - lNorm*( yc(nlat+1:2*nlat-2) - 2._dl*yc(nlat+2:2*nlat-1) + yc(nlat+3:2*nlat) )
    yp(3*nlat+2:4*nlat-1) = yp(3*nlat+2:4*nlat-1) + lNorm*( yc(2*nlat+1:3*nlat-2) - 2._dl*yc(2*nlat+2:3*nlat-1) + yc(2*nlat+3:3*nlat) )
    yp(2*nlat+2:3*nlat-1) = yp(2*nlat+2:3*nlat-1) - lNorm*( yc(3*nlat+1:4*nlat-2) - 2._dl*yc(3*nlat+2:4*nlat-1) + yc(3*nlat+3:4*nlat) )
#endif
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(R1)
    call laplacian_1d_wtype(tPair,dk)
    yp(I1) = yp(I1) + 0.5_dl*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I1)
    call laplacian_1d_wtype(tPair,dk)
    yp(R1) = yp(R1) - 0.5_dl*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(R2)
    call laplacian_1d_wtype(tPair,dk)
    yp(I2) = yp(I2) + 0.5_dl*tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I2)
    call laplacian_1d_wtype(tPair,dk)
    yp(R2) = yp(R2) - 0.5_dl*tPair%realSpace(:)
#endif
#endif
  end subroutine derivs

  !>@brief
  !> Solve coupled GPE for the simplified case of a symmetric 2 BEC experiment
  !>
  !> Solve the equations of motion for two symmetric coupled BECs.
  !> The Gross-Pitaevskii equations are solved in the linear complex basis, with the equations given by
  !> \f[
  !>    i\hbar\frac{d\psi_i}{dt} = \left-\frac{\hbar^2}{2m}\nabla^2 - g|\psi_i|^2 - g_c|\psi_{3-i}|^2 - \mu\right)\psi_i - \nu\psi_{3-j}
  !> \f]
  !> with the masses \f$m\f$ and self-couplings \f$g_{11} = g_{22} \equiv g\f$ assumed to be equal
  subroutine derivs_symmetric(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp    
  end subroutine derivs_symmetric
  
  subroutine derivs_angle_phase(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

  end subroutine derivs_angle_phase
  
  ! Equations of motion using Drummond's scaling of the variables
  subroutine derivs_drummond(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    real(dl), dimension(1:nLat,1:nFld) :: ysq

!    real(dl) :: nu_scl, rho_scl, del_m
    real(dl) :: c1, c2, c3

    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2
    
    c1 = sqrt(nu_scl)
    c2 = 0.5_dl*(nu_scl-1._dl)/c1
    c3 = 0.5_dl/(sqrt(nu_scl)*rho_scl)
    
    yp(R1) = 0._dl
    yp(I1) = 0._dl
    yp(R2) = 0._dl
    yp(I2) = 0._dl

    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    
    ! Fourier based differentiation
    ! Need to add Chebyshev in here
#ifdef DIFF
#ifdef FOURIER
    tPair%realSpace(:) = yc(R1)
    call laplacian_1d_wtype(tPair,dk)
    yp(I1) = yp(I1) + tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I1)
    call laplacian_1d_wtype(tPair,dk)
    yp(R1) = yp(R1) - tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(R2)
    call laplacian_1d_wtype(tPair,dk)
    yp(I2) = yp(I2) + tPair%realSpace(:)
    
    tPair%realSpace(:) = yc(I2)
    call laplacian_1d_wtype(tPair,dk)
    yp(R2) = yp(R2) - tPair%realSpace(:)
#endif
#endif
  end subroutine derivs_drummond
  
end module eom
