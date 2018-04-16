!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the variables and equations describing the model we are solving
!>
!> This module provides storage and equations of motion for a relativistic scalar
!> field evolving in one-spatial dimension.
!> The effective potential of the scalar is derived from a coupled cold atom system described
!> by the Gross-Pitaevskii equation, which is reflected in the various variable definitions
!> appearing below.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  implicit none

  integer, parameter :: nFld = 1, nLat = 512
  integer, parameter :: nVar = 2*nFld*nLat+1
  real(dl), dimension(1:nVar), target :: yvec

  real(dl), parameter :: nu = 2.e-3
  real(dl), parameter :: omega = 50._dl*2._dl*nu**0.5, del = (nu/2._dl)**0.5*(1._dl+0.5_dl)
  real(dl), parameter :: len = 0.5*50._dl / (2.*nu)**0.5, dx = len/dble(nLat), dk = twopi/len
  real(dl), parameter :: rho = 2.**(-3)*200._dl*2._dl*(nu)**0.5

  real(dl), parameter :: lambda = del*(2._dl/nu)**0.5
  real(dl), parameter :: m2eff = 4._dl*nu*(-1._dl+lambda**2)
  real(dl), parameter :: phi0 = 5._dl

#ifdef FOURIER
  type(transformPair1D) :: tPair
#endif
  
contains
  
  !>@brief
  !> Compute the time derivatives of the scalar with time-dependent coupling
  subroutine derivs_tdep(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
    real(dl) :: nueff
#ifdef DISCRETE
    real(dl), parameter :: lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    nueff = nu + del*omega*cos(omega*yc(TIND))

    yp(FLD) = yc(DFLD)
    yp(DFLD) = -4._dl*nueff*sin(yc(FLD))
    
#ifdef DIFF
    ! Discrete Differentiation
#ifdef DISCRETE
    ! Compute the endpoints
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#endif
    ! Fourier differentiation
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs_tdep

  !>@brief
  !> Compute the derivatives of the scalar field in the effective time-independent potential
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
    real(dl) :: lambda
#ifdef DISCRETE
    real(dl), parameter :: lNorm = 1._dl/dx**2
#endif
    lambda = del*(2./nu)**0.5
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -4._dl*nu*( sin(yc(FLD)) + 0.5_dl*lambda**2*sin(2._dl*yc(FLD)) )
    
#ifdef DIFF
    ! Discrete Differentiation
#ifdef DISCRETE
    ! Compute the endpoints
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#endif
    ! Fourier differentiation
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs

  
  !>@brief
  !> Compute the derivatives of the scalar field in the quadratically truncated potential
  subroutine derivs_free(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
    real(dl) :: lambda
#ifdef DISCRETE
    real(dl), parameter :: lNorm = 1._dl/dx**2
#endif
    lambda = del*(2./nu)**0.5
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -4._dl*nu*( -1._dl + lambda**2 )*yc(FLD)
    
#ifdef DIFF
    ! Discrete Differentiation
#ifdef DISCRETE
    ! Compute the endpoints
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#endif
    ! Fourier differentiation
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs_free

  !>@brief
  !> Equations of motion including the quantum pressure correction to the dispersion relation
  subroutine derivs_qp(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
    real(dl) :: nueff
    
#ifdef DISCRETE
    real(dl), parameter :: lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    nueff = nu + del*omega*cos(omega*yc(TIND))

    yp(FLD) = yc(DFLD)
    yp(DFLD) = -4._dl*nueff*sin(yc(FLD))
    
#ifdef DIFF
#ifdef DISCRETE
    ! Compute the endpoints
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#endif
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
    ! call laplacian_1d_wtype(tPair,dk)
    ! yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs_qp
  
  
end module eom
