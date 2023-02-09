!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>@author
!> Jonathan Braden
!> University College London
!>
!>@brief
!> Store the variables and equations describing the GPE we are solving
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
!#define DISCRETE T
#define FOURIER_DIFF T
#define DIFF T

#define R1 1:nLat
#define I1 (nLat+1):(2*nLat)
#define R2 (2*nLat+1):(3*nLat)
#define I2 (3*nLat+1):(4*nLat)
#define TIND 4*nLat+1

module eom
  use constants
#ifdef FOURIER
  use fftw3
#endif
  use gpe_model
  implicit none
    
contains
  
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
  
#ifdef ADJUST
  !>@brief
  !> Set the parameters of the Bose-Condensates assuming the simplified situation of equal couplings
  subroutine set_model_params_simple(mu_i,g_i,gc_i,nu_i)
    real(dl), intent(in), optional :: mu_i,g_i,gc_i,nu_i
    integer :: i
    
    if ( present(mu_i) ) mu = mu_i
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
  
  !>@brief
  !> Compute the derivatives of the vector of fields
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
    real(dl), dimension(1:nLat,1:nFld) :: ysq  ! probably faster to preallocate this
#ifdef DISCRETE
    real(dl) :: lNorm
    lNorm = 0.5_dl/dx**2
#endif
    
    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2
    yp(TIND) = 1._dl
    
    nueff = nu + del*omega*cos(omega*yc(TIND))
    
    yp(R1) = -mu*yc(I1) + ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(I1)
    yp(I1) = mu*yc(R1)  - ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(R1)
    yp(R2) = -mu*yc(I2) + ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(I2)
    yp(I2) = mu*yc(R2)  - ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(R2)

    yp(R1) = yp(R1) - nueff*yc(I2)
    yp(I1) = yp(I1) + nueff*yc(R2)
    yp(R2) = yp(R2) - nueff*yc(I1)
    yp(I2) = yp(I2) + nueff*yc(R1)
    
#ifdef DIFF
#ifdef DISCRETE
!!! Discrete Derivatives
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
!!! Fourier Derivatives
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
#ifdef CHEBY_DIFF
    
#endif
#endif
  end subroutine derivs

  subroutine derivs_no_oscillation(yc, yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    real(dl), dimension(1:nLat,1:nFld) :: ysq
    
    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2
    yp(TIND) = 1._dl

!    yp(R1) = -mu*yc(I1) + ()*yc(I1)
!    yp(I1) = mu*yc(R1) - ()*yc(R1)
!    yp(R2) = -mu*yc(I2) + ()*yc(I2)
!    yp(I2) = mu*yc(R2) - ( )*yc(R2)

  end subroutine derivs_no_oscillation
  
  !>@brief
  !> Solve coupled GPE for the simplified case of a symmetric 2 BEC experiment
  !>
  !> Solve the equations of motion for two symmetric coupled BECs.
  !> The Gross-Pitaevskii equations are solved in the linear complex basis, with the equations given by
  !> \f[
  !>    i\hbar\frac{d\psi_i}{dt} = \left-\frac{\hbar^2}{2m}\nabla^2 - g|\psi_i|^2 - g_c|\psi_{3-i}|^2 - \mu\right)\psi_i - \nu\psi_{3-j}
  !> \f]
  !> with the masses \f$m\f$ and self-couplings \f$g_{11} = g_{22} \equiv g\f$ assumed to be equal
  subroutine derivs_0d(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    real(dl) :: nueff
    real(dl), dimension(1:nLat,1:nFld) :: ysq
 
    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2

    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    nueff = nu + del*omega*cos(omega*yc(TIND))
    
    ! These are ugly, add some default vectors so vectorisation can be done more easily
    yp(R1) = -mu*yc(I1) + ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(I1) 
    yp(I1) = mu*yc(R1)  - ( lam(1,1)*ysq(:,1) + lam(1,2)*ysq(:,2) )*yc(R1) 
    yp(R2) = -mu*yc(I2) + ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(I2) 
    yp(I2) = mu*yc(R2)  - ( lam(2,1)*ysq(:,1) + lam(2,2)*ysq(:,2) )*yc(R2)

    yp(R1) = yp(R1) - nueff*yc(I2)
    yp(I1) = yp(I1) + nueff*yc(R2)
    yp(R2) = yp(R2) - nueff*yc(I1)
    yp(I2) = yp(I2) + nueff*yc(R1)
  end subroutine derivs_0d
  
  subroutine derivs_0d_symmetric(yc,yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp

    real(dl) :: nueff
    real(dl), dimension(1:nLat,1:nFld) :: ysq

    
    ysq(:,1) = yc(R1)**2 + yc(I1)**2
    ysq(:,2) = yc(R2)**2 + yc(I2)**2

    nueff = nu + del*omega*cos(omega*yc(TIND))

    yp(R1) = -mu*yc(I1) + gs*ysq(:,1)*yc(I1) - nueff*yc(I2)
    yp(I1) = mu*yc(R1)  - gs*ysq(:,1)*yc(R1) + nueff*yc(R2)
    yp(R2) = -mu*yc(I2) + gs*ysq(:,2)*yc(I2) - nueff*yc(I1)
    yp(I2) = mu*yc(R2)  - gs*ysq(:,2)*yc(R2) + nueff*yc(R1)

    yp(TIND) = 1._dl
  end subroutine derivs_0d_symmetric
    
  !>@brief
  !> Compute the energy density of the underlying condensates, using the correctly defined gradient energy
  subroutine compute_energy(rho,fld)
    real(dl), dimension(1:nlat), intent(out) :: rho
    real(dl), dimension(1:nFld,1:nLat), intent(in) :: fld
    rho = 0._dl
#ifdef DIFF
#ifdef DISCRETE
    rho = rho + 0._dl
    rho = rho + 0._dl
    rho = rho + 0._dl
    rho = rho + 0._dl
#endif
#ifdef FOURIER_DIFF
    tPair%realSpace(:) = fld(1,:); call gradsquared_1d_wtype(tPair,dk)
    rho = rho + tPair%realSpace(:)
    tPair%realSpace(:) = fld(2,:); call gradsquared_1d_wtype(tPair,dk)
    rho = rho + tPair%realSpace(:)
    tPair%realSpace(:) = fld(3,:); call gradsquared_1d_wtype(tPair,dk)
    rho = rho + tPair%realSpace(:)
    tPair%realSpace(:) = fld(4,:); call gradsquared_1d_wtype(tPair,dk)
    rho = rho + tPair%realSpace(:)
#endif
#endif
    !rho = rho + 0.5_dl*lam(1,1)*(yc(R1)**1+yc(I1)**2)**2 + 0.5_dl*lam(2,2)*(yc(R2)**2+yc(I2)**2)**2 + lam(1,2)*(yc(R1)**2+yc(I1)**2)*(yc(R2)**2+yc(I2)**2) - 2._dl*nu*(yc(R1)*yc(R2)+yc(I1)*yc(I2))
  end subroutine compute_energy

  ! This needs to be optimized
  real(dl) function compute_total_energy(fld) result(en)
    real(dl), dimension(:,:,:), intent(in) :: fld

    integer :: nfld
    integer :: i,j

    nfld = size(fld,dim=3)
    en = 0._dl
#ifdef DIFF
#ifdef DISCRETE
#endif
#ifdef FOURIER_DIFF
    do i=1,nfld
       do j=1,2  ! Real and complex components
          tPair%realSpace(:) = fld(:,j,i)
          call laplacian_1d_wtype(tPair,dk)
          en = en - 0.5_dl*sum(fld(:,j,i)*tPair%realSpace)
          !call gradsquared_1d_wtype(tPair,dk)
          !en = en + 0.5_dl*sum(tPair%realSpace)
       enddo
    enddo
#endif
#endif
    do i=1,nfld
       do j=1,nfld
          en = en + 0.5_dl*lam(i,j)*sum( (fld(:,1,i)**2+fld(:,2,i)**2)*(fld(:,1,j)**2+fld(:,2,j)**2) )
       enddo
    enddo

    do i=1,nfld
       do j=1,nfld
          en = en - nu*sum( fld(:,1,i)*fld(:,1,j) + fld(:,2,i)*fld(:,2,j) )
       enddo
    enddo

    en = en / size(fld,dim=1)
  end function compute_total_energy

  real(dl) function compute_grad_energy(fld) result(ge)
    real(dl), dimension(:,:,:), intent(in) :: fld

    integer :: i,j

    ge = 0._dl
    do i=1,2
       do j=1,2
          tPair%realSpace(:) = fld(:,j,i)
          call gradsquared_1d_wtype(tPair,dk)
          ge = ge + 0.5_dl*sum(tPair%realSpace)
       enddo
    enddo
  end function compute_grad_energy

  real(dl) function compute_scattering_energy(fld) result(pe)
    real(dl), dimension(:,:,:), intent(in) :: fld

    integer :: nfld, i, j

    nfld = size(fld,dim=3)
    pe = 0._dl
    do i=1,nfld
       do j=1,nfld
          pe = pe + 0.5_dl*lam(i,j)*sum( (fld(:,1,i)**2+fld(:,2,i)**2)*fld(:,1,j)**2+fld(:,2,j)**2 )
       enddo
    enddo
  end function compute_scattering_energy
  
  real(dl) function compute_total_density(fld) result(rho)
    real(dl), dimension(:,:,:), intent(in) :: fld
    rho = sum(fld**2) / size(fld,dim=1)
  end function compute_total_density
  
end module eom
