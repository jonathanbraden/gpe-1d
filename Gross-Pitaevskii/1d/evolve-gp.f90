program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use gaussianRandomField
  use eom
  use integrator
  implicit none
  real(dl), dimension(:,:,:), pointer :: fld
  
  fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  call initialise_fields(fld)
  call setup(nVar)
  call output_fields(fld)
!  call time_evolve(0.2_dl/omega,80000,800)
  call time_evolve(0.2_dl/omega,3200,800)
  
contains

  !>@brief
  !> Compute the desired time step dt adaptively using precomputed conditions.
  !>
  !>@todo
  !>@arg Write this.
  real(dl) function get_dt() result(dt)

    ! Need to satisfy:
    ! - Courant Condition
    ! - Resolve short time oscillator
    ! - Characteristic oscillation frequency
  end function get_dt
  
  !>@brief
  !> Initialise the integrator, setup FFTW, boot MPI, and perform other necessary setup
  !> before starting the program
  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
    call initialize_transform_1d(tPair,nLat)
  end subroutine setup

  !>@brief
  !> Find the false vacuum minima around which we will evolve our field fluctuations initially
  subroutine initialise_mean_fields(fld,rho)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), intent(in) :: rho
    fld(:,1,1) =  rho; fld(:,2,1) = 0._dl
    fld(:,1,2) = -rho; fld(:,2,2) = 0._dl
  end subroutine initialise_mean_fields

  subroutine initialise_fluctuations(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl) :: df(1:nLat), spec(1:nLat/2+1)
    integer :: i,j
    
    spec = 0._dl
    !    do i=2,nLat/4
    do i=nLat/4-4,nLat/4
       spec(i) = 100._dl / (sqrt(len*rho))
    enddo
    do i = 1,nFld; do j=1,2
       call generate_1dGRF(df,spec,.false.)
       fld(:,i,j) = fld(:,i,j) + df
    enddo; enddo
  end subroutine initialise_fluctuations
  
  subroutine initialise_fields(fld)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), parameter :: rho = 1._dl
    integer :: i; real(dl) :: dt, theta
    
    call initialise_mean_fields(fld,rho)
    yvec(4*nLat+1) = 0._dl ! Add a tcur pointer here
!#define SW_INIT T
#ifdef SW_INIT
    dt = twopi / dble(nLat)
    fld(:,1,1) = rho; fld(:,2,1) = 0._dl
    do i=1,nLat
       theta = (i-1)*dt
       fld(i,1,2) = rho*cos(theta); fld(i,2,2) = rho*sin(theta)
    enddo
#else
    call initialise_fluctuations(fld)
#endif
  end subroutine initialise_fields
  
  subroutine time_evolve(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no
    integer :: i,j, outsize, nums

    print*,"dx is ", dx, "dt is ",dt
    if (dt > dx**2) print*,"Warning, violating Courant condition"
    
    outsize = ns/no; nums = ns/outsize
    print*,"dt out is ",dt*outsize
    do i=1,nums
       do j=1,outsize
          call gl10(yvec,dt)
       enddo
       call output_fields(fld)
    enddo
  end subroutine time_evolve

  subroutine output_fields(fld)
    real(dl), dimension(:,:,:), intent(in) :: fld
    logical :: o; integer :: i

    inquire(file='fields.dat',opened=o)
    if (.not.o) open(unit=99,file='fields.dat')

    ! Fix this if I change array orderings
    do i=1,nLat
       write(99,*) fld(i,:,:)
    enddo
    write(99,*)
    
  end subroutine output_fields
  
end program Gross_Pitaevskii_1d
