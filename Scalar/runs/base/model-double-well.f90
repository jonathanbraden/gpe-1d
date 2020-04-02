module model
  use constants, only : dl, twopi
  implicit none
  real(dl) :: eps, phi0

contains
  subroutine set_model(eps_,phi0_)
    real(dl), intent(in) :: eps_, phi0_
    eps = eps_; phi0 = phi0_
  end subroutine set_model

  elemental function v(phi)
    real(dl) :: v
    real(dl), intent(in) :: phi
    v = 0.25_dl*(phi**2-1._dl)**2 + eps*(phi/3._dl-1._dl)
  end subroutine v

  elemental function vp(phi)
    real(dl) :: vp
    real(dl), intent(in) :: phi
    vp = (phi**2-1._dl)*(phi+eps)
  end function vp

  elemental function vpp(phi)
    real(dl) :: vpp
    real(dl), intent(in) :: phi
    vpp = 3._dl*phi**2-1._dl + 2._dl*eps*phi
  end function vpp

  real(dl) function phi_fv()
    phi_fv = -1._dl
  end function phi_fv

  real(dl) function phi_tv()
    phi_tv = 1._dl
  end function phi_tv

end module model
