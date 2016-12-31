program Gross_Pitaevskii_1d
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use eom
  use integrator
  implicit none
  real(dl), dimension(:,:,:), pointer :: fld
  
  fld(1:nLat,1:2,1:nFld) => yvec(1:2*nLat*nFld)
  call initialise_fields(fld)
  call setup(nVar)
  call output_fields(fld)
  call time_evolve(0.05_dl,10000,1000)
  
contains

  !>@brief
  !> Compute the desired time step dt adaptively using precomputed conditions.
  !>
  !>@todo
  !>@arg Write this.
  real(dl) function get_dt() result(dt)
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
  subroutine initialise_mean_fields()
  end subroutine initialise_mean_fields
  
  subroutine initialise_fields(fld)
    real(dl), dimension(:,:,:), intent(out) :: fld
    real(dl), parameter :: rho = 1._dl
    integer :: i; real(dl) :: dt, theta
    
    call initialise_mean_fields()
    dt = twopi / dble(nLat)
    fld(:,1,1) = rho; fld(:,2,1) = 0._dl
    do i=1,nLat
       theta = (i-1)*dt
       fld(i,1,2) = rho*cos(theta); fld(i,2,2) = rho*sin(theta)
    enddo
    !call initialise_fluctuations()
  end subroutine initialise_fields
  
  subroutine time_evolve(dt,ns,no)
    real(dl), intent(in) :: dt
    integer, intent(in) :: ns, no
    integer :: i,j, outsize, nums

    print*,"dx is ", dx
    if (dt > dx**2) print*,"Warning, violating Courant condition"
    
    outsize = ns/no; nums = ns/outsize
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
