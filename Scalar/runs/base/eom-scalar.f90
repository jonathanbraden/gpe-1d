!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> This module provides storage and equations of motion for a relativistic scalar evolving in one spatial dimension
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

  ! Fix these so that they aren't fixed.
!  integer, parameter :: nLat=512, nFld=1
!  integer, parameter :: nVar = 2*nLat*nFld+1
!  real(dl), dimension(1:nVar), target :: yvec  ! This is a problem, does it work if I make it allocatable?

  integer :: nLat, nFld, nVar
  real(dl), dimension(:), allocatable, target :: yvec

  real(dl) :: len, dx, dk
  real(dl) :: lambda, m2eff

#ifdef FOURIER
  type(transformPair1D) :: tPair
#endif
  
contains

  ! Fix the paramter nature above (in particular the declaration of yvec
  subroutine set_lattice_params(n,l,nf)
    integer, intent(in) :: n, nf
    real(dl), intent(in) :: l
    nLat = n; len = l; nFld = nf

    nVar = 2*nFld*nLat + 1; allocate(yvec(1:nVar))
    dx = len / dble(nLat); dk = twopi/len
  end subroutine set_lattice_params

  ! Add appropriate subroutine calls here to potential derivs, etc. here
  subroutine set_model_params(m2,lam)
    real(dl), intent(in) :: m2, lam
    lambda = lam; m2eff = m2*(-1._dl+lambda**2)
  end subroutine set_model_params

  !>@brief
  !> Compute the derivatives of the scalar field in the effective time-independent potential
  subroutine derivs(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
#ifdef DISCRETE
    real(dl), parameter :: lNorm = 1._dl/dx**2
#endif
    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)
    yp(DFLD) = -( sin(yc(FLD)) + 0.5_dl*lambda**2*sin(2._dl*yc(FLD)) )
    
#ifdef DIFF
#ifdef DISCRETE
    yp(nlat+1) = yp(nlat+1) + lNorm * ( yc(nlat) - 2._dl*yc(1) + yc(2) )
    yp(2*nlat) = yp(2*nlat) + lNorm * ( yc(nlat-1) - 2._dl*yc(nlat) + yc(1) )
    
    yp(nlat+2:2*nlat-1) = yp(nlat+2:2*nlat-1) + lNorm*( yc(1:nlat-2) - 2._dl*yc(2:nlat-1) + yc(3:nlat) )
#endif
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
#endif
#endif
  end subroutine derivs
  
end module eom
