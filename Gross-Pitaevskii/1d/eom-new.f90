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

#define FOURIER T
!#define DIFF T
#define FIELD_TYPE Field_Model
module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  implicit none

  integer, parameter :: nFld = 2, nLat = 32
  integer, parameter :: nVar = 2*nFld*nLat
  real(dl), dimension(1:nVar), target :: yvec
  real(dl), parameter :: len = 80._dl, dx = len/dble(nLat), dk = twopi/len
  real(dl), parameter :: mu = 0.1_dl, gs = 0.1_dl, gc = 0._dl, nu = 0.1_dl
  real(dl), dimension(1:nFld,1:nFld), parameter :: lam = reshape( (/ gs, gc, gc, gs /), [nFld,nFld] )

#ifdef FOURIER
  type(transformPair1D) :: tPair
#endif

  type Grid_Params
     integer :: nLat
     real(dl) :: len, dx, dk
  end type Grid_Params

  type Model_Params
     integer :: nFld
     real(dl), dimension(:,:), allocatable :: lam, nu
     real(dl), dimension(:), allocatable :: mu, mass
  end type Model_Params
  
  type FIELD_TYPE
     integer :: nx
     integer :: nfld, nvar
     real(dl) :: len, dx, dk
     real(dl), dimension(:), allocatable :: yvec
     ! Add a pointer to reshape the variables for ease of coding
  end type FIELD_TYPE

  type(Model_Params) :: params
  
contains

#ifdef ADJUST
  !>@brief
  !> Set the parameters of the Bose-Condensates assuming the simplified situation of equal couplings
  subroutine set_model_params_simple(mu_i,g_i,gc_i,nu_i)
    real(dl), intent(in), optional :: mu_i,g_i,gc_i,nu_i
    integer :: i
    
    if ( present(mu_i) )       mu = mu_i
    if ( present(g_i) )  then; do i=1,nFld; lam(i,i) = g_i; enddo; endif
    if ( present(gc_i) ) then; endif  ! Finish this one
    if ( present(nu_i) ) then; nu = nu_i; do i=1,nFld; nu(i,i) = 0.1_dl; enddo; endif
  end subroutine set_model_params_simple
  
  subroutine set_model_params(mu_i,g_i,nu_i)
    real(dl), optional, intent(in) :: mu_i, g_i(1:nFld,1:nFld), nu_i(1:nFld,1:nFld)

    if (present(mu_i)) mu = mu_i
    if (present(g_i)) lam = g_i
    if (present(nu_i)) nu = nu_i
  end subroutine set_model_params
#endif
  
  type(Model_Params) function new_model_params(nf) result(this)
    integer, intent(in) :: nf
    allocate( this%lam(1:nf,1:nf),this%nu(1:nf,1:nf) )
    allocate( this%mu(1:nf), this%mass(1:nf) )
  end function new_model_params
  
  subroutine create_model_params(this,nf)
    type(Model_Params), intent(out) :: this
    integer, intent(in) :: nf
    this%nFld = nf
    allocate( this%lam(1:nf,1:nf), this%nu(1:nf,1:nf) )
    allocate( this%mu(1:nf) )
  end subroutine create_model_params
  
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

  ! Rather than hardcoding parameters above, use this subroutine
  subroutine initialise_field(n)
    integer, intent(in) :: n
    !call initialise_transform_1d(tPair,n)
  end subroutine initialise_field
  
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
    
    real(dl), dimension(1:nLat,1:nFld) :: ysq
    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2

!    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    
    ! These are ugly, add some default vectors so vectorisation can be done more easily
    yp(R1) = -mu*yc(I1) + ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(I1) ! -laplacian(I1)
    yp(I1) = mu*yc(R1)  - ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(R1) ! +laplacian(R1)
    yp(R2) = -mu*yc(I2) + ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(I2)  ! -laplacian(I2)
    yp(I2) = mu*yc(R2)  - ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(R2) ! + laplacian(R2)

    yp(R1) = yp(R1) - nu*yc(I2)
    yp(I1) = yp(I1) + nu*yc(R2)
    yp(R2) = yp(R2) - nu*yc(I1)
    yp(I2) = yp(I2) + nu*yc(R1)

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
  end subroutine derivs
  
end module eom
