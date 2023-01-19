module HomogeneousField
  use constants, only : dl, twopi
  
contains

  !>@brief
  !> Induce a global phase rotation of a field.
  !
  !>@param[inout] fld : The fields with shape (nLat, 1:2, nFld)
  !>@param[in] phi : The relative phase to imprint
  !>@param[in] i : Index of field to rotate
  subroutine rotate_phase(fld, phi, i)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: phi
    integer, intent(in) :: i

    real(dl) :: fld_tmp(1:size(fld,dim=1))
    real(dl) :: cphi, sphi

    cphi = cos(phi); sphi = sin(phi)
    fld_tmp = fld(:,1,i)
    
    fld(:,1,i) = fld(:,1,i)*cphi - fld(:,2,i)*sphi
    fld(:,2,i) = fld_tmp*sphi + fld(:,2,i)*cphi
  end subroutine rotate_phase

  ! This can be optimized to not call another subroutine
  subroutine rotate_global_phase(fld, theta)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: theta

    integer :: i, num_fld

    num_fld = size(fld,dim=3)
    do i=1,num_fld
       call rotate_phase(fld, theta, i)
    enddo
  end subroutine rotate_global_phase

  subroutine imprint_relative_phase(fld, dphi, i, j)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: dphi
    integer, intent(in) :: i,j

    call rotate_phase(fld, 0.5_dl*dphi, j)
    call rotate_phase(fld, -0.5_dl*dphi, i)
  end subroutine imprint_relative_phase

  subroutine shift_homogeneous_field(fld, shift)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), dimension(1:2,1:size(fld,dim=3)), intent(in) :: shift

    integer :: i,j

    do i=1,size(fld,dim=3)
       do j=1,2
          fld(:,j,i) = fld(:,j,i) - shift(j,i)
       enddo
    enddo
  end subroutine shift_homogeneous_field
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solution of homogeneous FV state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    real(dl) function find_max_angle(nu,dg,gc) result(theta)
    real(dl), intent(in) :: nu, dg, gc
    real(dl) :: nu_bg, dg_bg, dth; integer :: l
    integer, parameter :: maxit = 16; real(dl), parameter :: eps = 1.e-16

    theta = 0.125_dl*twopi + 0.25_dl*dg/(1._dl-gc - nu)
    nu_bg = nu / (1._dl-gc); dg_bg = dg / (1._dl-gc)

    do l = 1,maxit
       dth = - dvTheta(theta,nu_bg,dg_bg) / d2vTheta(theta,nu_bg,dg_bg)
       theta = theta + dth
       if (abs(dth) < eps) exit
    end do
    
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Convenience functions to compute things like dispersion relations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fix to allow for non-unit g values, etc
  real(dl) function dispersion_total_symmetric(k) result(om)
    real(dl), intent(in) :: k

    om = k*sqrt(1._dl + 0.25_dl*k**2)
  end function dispersion_total_symmetric

  ! Fix this to allow for non-unit g values
  real(dl) function dispersion_relative_symmetric(k,nu) result(om)
    real(dl), intent(in) :: k, nu

    om = 2._dl*sqrt( (nu + 0.25_dl*k**2) * (1._dl + nu + 0.25_dl*k**2) )
  end function dispersion_relative_symmetric
  
end module HomogeneousField
