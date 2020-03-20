module fluctuations
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use gaussianRandomField
  
  implicit none
  complex(C_DOUBLE_COMPLEX), parameter :: iImag = (0._dl,1._dl)  ! define this as a private constant

  !>@brief
  !> Basic interface for writing function to return a spectrum
  abstract interface
     real(dl) function spec_null(k)
       use constants, only : dl
       real(dl), intent(in) :: k
     end function spec_null
  end interface

  !>@brief
  !> Basic interface for writing full correlation information
  abstract interface
     function corr_null(k,nf)
       use constants, only : dl
       real(dl), dimension(1:nf,1:nf) :: corr_fun
       real(dl), intent(in) :: k
       integer, intent(in) :: nf
     end function corr_null
  end interface

contains

  !>@brief
  !> Generate a realisation of a GRF with the spectrum specified by the 
  !> input function
  subroutine initial_fluctuations(dfld,spec)
    real(dl), dimension(:,:), intent(out) :: dfld
    procedure(spec_null) :: spec
  end subroutine initial_fluctuations

  !>@brief
  !> Generate a set of correlated GRFs with specified covariance matrix
  !> in spectral space.
  subroutine initial_fluctuations_correlated(dfld,corr,nf,n,dk)
    real(dl), dimension(:,:), intent(out) :: dfld
    integer, intent(in) :: nf, n
    real(dl), intent(in) :: dk
    procedure(corr_null) :: corr

    integer :: i, nn
    real(dl) :: norm, k
    complex(C_DOUBLE_COMPLEX) :: Fk(1:n/2+1,1:nf), deviate
    real(dl), dimension(1:nf,1:nf) :: chol
    real(dl), dimension(1:nf) :: amp, phase

    ! Need to initialise random numbers here

    nn = n/2+1
    norm = 1._dl  ! Fix this
    do i=1,nn
       k = dk*(i-1); chol = corr(k,nf)
       ! call DPOTRF('L')   ! Do Cholesky decomposition
       call random_number(amp(1:nf)); call random_number(phase(1:nf))
       Fk(i,:) = matmul(chol,sqrt(-log(amp))*exp(iImag*twopi*phase)
    enddo

    ! Fourier transform the fields back.  Worth doing multiple in parallel?
  end subroutine initial_fluctuations_correlated

  !>@brief
  !> Create a spectrum suitable for use in GRF generation from a function specifying the spectrum
  function create_spectrum(spec,n,dk) result(s)
    real(dl), dimension(1:n/2+1) :: s
    procedure(spec_null) :: spec
    integer, intent(in) :: n
    real(dl), intent(in) :: dk

    integer :: i
    do i=1,n/2+1; s(i) = spec((i-1)*dk); enddo
  end function create_spectrum

  subroutine initialise_fluctuations(fld,type,kmax,phi0)
    real(dl), dimension(:,:), intent(inout) :: fld
    integer, intent(in) :: type  ! Turn this into a string
    integer, intent(in), optional :: kmax

    integer :: km

    ! Need to include spec
    km = size(spec); if (present(kmax)) km = kmax
    phiL = 0.5_dl; if (present(phi0)) phiL = phi0

    select case(type)
       case(1)
          ! Minkowski
       case(2)
          ! Thermal
       case(3)
          ! Thermal High-T
  end subroutine initialise_fluctuations

  subroutine initialise_fluctuations_vacuum(fld)
    real(dl), dimension(:,:), intent(inout) :: fld
  end subroutine initialise_fluctuations_vacuum

  subroutine initialise_fluctuations_thermal(fld)
    real(dl), dimension(:,:), intent(inout) :: fld
  end subroutine initialise_fluctuations_thermal

  subroutine initialise_fluctuations_high_temp()
  end subroutine initialise_fluctuations_high_temp

  subroutine initialise_fluctuations_fields()
  end subroutine initialise_fluctuations_fields

  subroutine initialise_fluctuations_momenta()
  end subroutine initialise_fluctuations_momentua

! Temporary GRF generation goes here.  Once debugged, move into separate module
  subroutine correlated_GRF_1D_real(df,nlat,nf,spec)
    real(dl), dimension(1:n,1:nf), intent(out) :: df
    integer, intent(in) :: nlat, nf
    real(dl), intent(in), dimension(1:nf,1:nf) :: spec  ! Make this a function I pass in
    
    integer :: nn, nnk
    real(dl) :: amp, phase
    complex(dl) :: deviate
    complex(C_DOUBLE_COMPLEX), allocatable :: Fk(:,:,:)
    type(C_PTR) :: fft_plan

    integer :: i
    integer :: info
    ! Storage for CHOLESKY decompostion
    real(dl), dimension(1:nf,1:nf) :: chol

    call initialize_rand(75,13)  !!! To Do: Actually make this random

    nn = nlat/2+1
    nnk = nn  ! Modify to allow for spectral cuts
    allocate(Fk(1:nn,1:nf,1:nf))

    ! Do either Cholesky decomposition or eigendecomposition of correlation function
    ! DPOTRF is LAPACK routine
    ! DPPTRF used packed storage

    ! As a check, can also do diagonalisation of matrix

    do i=2,nnk
       chol = spec(i)
       ! call DPOTRF('L',nf,chol,nf,info)  ! Check upper vs lower
       ! if (info .ne. 0) ! Check for errors
       call random_number(amp); call random_number(phase)
       ! Get complex Gaussian realisations via random numbers
       deviate = sqrt(-log(amp))*exp(iImag*twopi*phase)     ! Check normalisation here
       Fk(i,1:nf) = matmul(chol,deviate)
    enddo
    Fk(1) = 0._dl
    Fk(nnk+1:nn) = 0._dl

    fft_plan = fftw_plan_dft_c2r_1d(nlat,Fk,df,FFTW_ESTIMATE)  ! This is wrong for the multifield FFT
    call fftw_execute_dft_c2r(fft_plan,Fk,df)                  ! Also wrong for multifield FFT
  end subroutine correlated_GRF_1D

  subroutine convolve_GRF_1D(df,nlat,nf)
    real(dl), dimension(1:n,1:nf), intent(out) :: df
    integer, intent(in) :: nlat, nf
  end subroutine convolve_GRF_1D

end module fluctuations
