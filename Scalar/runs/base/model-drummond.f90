module model
  use constants, only : dl, twopi
  implicit none
  real(dl) :: lambda, phi0

contains
  subroutine set_model(lam_,phi0_)
    real(dl), intent(in) :: lam_, phi0_
    lambda = lam_; phi0 = phi0_
  end subroutine set_model

  elemental function v(phi)
    real(dl) :: v
    real(dl), intent(in) :: phi
    v = 1._dl - cos(phi) + 0.5_dl*lambda**2*sin(phi)**2
  end function v

  elemental function vp(phi)
    real(dl) :: vp
    real(dl), intent(in) :: phi
    vp = sin(phi) + 0.5_dl*lambda**2*sin(2._dl*phi)
  end function vp

  elemental function vpp(phi)
    real(dl) :: vpp
    real(dl), intent(in) :: phi
    vpp = cos(phi) + lambda**2*cos(2._dl*phi)
  end function vpp

  real(dl) function phi_fv()
    phi_fv = 0.5*twopi
  end function phi_fv

  real(dl) function phi_tv()
    phi_tv = 0._dl
  end function phi_tv

end module model
