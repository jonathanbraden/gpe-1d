! TO DO : Remove the len nonlocality (and any others) in this code
! TO DO : Finish the wrapper function

module Fluctuations
  use constants, only : dl, twopi
  use gaussianRandomField
  use gpe_model ! remove this dependence. It stores grid params, etc.
  
  implicit none

  type SpecParams
     real(dl) :: rho, len, m2eff  ! Get rid of len
     real(dl) :: nu, lamEff
     integer :: nCut = 1
     character(8) :: type = 'BOGO'
     logical, dimension(2) :: modes = (/.true.,.true./) ! Update for nFlds
  end type SpecParams
  
contains

  ! Get this to take in model parameters
  function make_spec_params(rho, len, nu, lamEff, type, modes, nCut) result(params)
    real(dl), intent(in) :: rho, len, nu, lamEff
    character(*), intent(in) :: type
    logical, dimension(1:2), intent(in) :: modes
    integer, intent(in) :: nCut
    
    type(SpecParams) :: params

    params%nCut = nCut
    params%rho = rho
    params%len = len
    params%nu = nu
    params%lamEff = lamEff
    params%m2eff = 4._dl*nu*(lamEff**2-1._dl)
    params%type = trim(type)
    params%modes(1:2) = modes(1:2)
  end function make_spec_params

  subroutine print_spec_params(params)
    type(SpecParams), intent(in) :: params

    print*,"===================================="
    print*,"Fluctuation Properties"
    print*,"------------------------------------"
    print*,"type of fluctuations   : ", params%type
    print*,"m^2_eff                : ", params%m2eff
    print*,"Initialise total modes : ", params%modes(1)
    print*,"Initialise rel modes   : ", params%modes(2)
    print*,"Spectral cutff number  : ", params%nCut
    print*,"===================================="
  end subroutine print_spec_params
  
  subroutine initialise_fluctuations(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    select case (params%type)
    case ('KG')
       call initialise_fluctuations_kg(fld, params)
    case ('BOGO')
       call initialise_fluctuations_bogoliubov(fld, params)
    case ('WHITE')
       call initialise_fluctuations_white(fld, params)
    case default
       call initialise_fluctuations_white(fld, params)
    end select
  end subroutine initialise_fluctuations

  subroutine initialise_fluctuations_white(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    real(dl), dimension( 1:size(fld,dim=1) ) :: df
    real(dl), dimension( 1:size(fld,dim=1)/2+1 ) :: spec
    integer :: num_fld, n
    integer :: i, j

    real(dl) :: rho, len
    integer :: nCut

    rho = params%rho
    len = params%len
    nCut = params%nCut
    
    n = size(fld,dim=1)
    num_fld = size(fld,dim=3)
    
    spec = 0._dl
    do i=2,n/2; do j=1,2
       spec(i) = 1._dl / sqrt(2._dl*len*rho) ! Change this version
    enddo; enddo

    do i = 1,num_fld; do j=1,2
       call generate_1dGRF(df, spec(1:nCut), .false.)
       fld(:,j,i) = fld(:,j,i) + df
    enddo; enddo
  end subroutine initialise_fluctuations_white

  subroutine initialise_fluctuations_kg(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams) :: params

    real(dl), dimension(1:size(fld,dim=1),1:2) :: df_rel, df_tot
    real(dl) :: spec(1:size(fld,dim=1)/2+1)
    real(dl) :: m2eff, dk, keff, norm
    integer :: nCut, i,j
    logical, dimension(2) :: modes
    real(dl) :: len, rho

    len = params%len
    rho = params%rho
    
    modes(1:2) = params%modes(1:2)    
    m2eff = params%m2eff
    nCut = params%nCut
    
    norm = 1._dl / sqrt(2._dl*len*rho)
    dk = twopi/len

    df_rel = 0._dl
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl/sqrt(keff)
    enddo
    spec = spec * norm
    do i = 1,nFld; do j=1,2
       if (modes(2)) call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo; enddo

    df_tot = 0.
    spec = 0._dl
    do i=1,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl / (keff**2+m2eff)**0.25 
    enddo
    spec = spec * norm
    do i=1,nFld; do j=1,2
       if (modes(1)) call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo; enddo

    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo
  end subroutine initialise_fluctuations_kg
  
  subroutine initialise_fluctuations_bogoliubov(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    real(dl), dimension(1:size(fld,dim=1),1:2) :: df_rel, df_tot
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec
    integer :: i,j
    real(dl) :: norm, dk, keff
    real(dl) :: lameff, nu_

    integer :: nCut

    df_rel = 0._dl
    df_tot = 0._dl
    
    nCut = params%nCut
    lameff = params%lamEff
    nu_ = params%nu
    
    ! Fix all this nonlocality
    norm = 1._dl / sqrt(2._dl*len*params%rho)
    dk = twopi / len   ! Fix this nonlocality
    
    ! Relative modes
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = sqrt( (keff**2+2._dl)/(keff*sqrt(keff**2+4._dl)) )
    enddo
    spec = spec * norm
    do j=1,2
       if (params%modes(2)) call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo
    
    ! Now get total modes
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = sqrt( (keff**2 + 2._dl - 4._dl*nu_)  / sqrt(keff**2+4._dl*nu_*(lameff**2-1._dl)) / sqrt(keff**2+4._dl-4._dl*nu_*(lameff**2+1._dl)) )
    enddo
    spec = spec * norm
    do j=1,2
       if (params%modes(1)) call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo
    
    ! Now remix the fields and add them to the background
    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo   
  end subroutine initialise_fluctuations_bogoliubov

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More verbose versions of calls, with parameters explicitly listed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine initialise_fluctuations_white_long(fld, rho, nCut)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: rho
    integer, intent(in) :: nCut

    real(dl), dimension( 1:size(fld,dim=1) ) :: df
    real(dl), dimension( 1:size(fld,dim=1)/2+1 ) :: spec
    integer :: num_fld, n
    integer :: i, j
    
    n = size(fld,dim=1)
    num_fld = size(fld,dim=3)
    
    spec = 0._dl
    do i=2,n/2; do j=1,2
       spec(i) = 1._dl / sqrt(2._dl*len*rho) ! Change nonlocality here
    enddo; enddo

    do i = 1,num_fld; do j=1,2
       call generate_1dGRF(df, spec(1:nCut), .false.)
       fld(:,j,i) = fld(:,j,i) + df
    enddo; enddo
  end subroutine initialise_fluctuations_white_long

  subroutine initialise_fluctuations_kg_long(fld, rho, m2eff, nCut, modes)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: rho, m2eff
    integer, intent(in) :: nCut
    logical, intent(in), dimension(1:size(fld,dim=3)) :: modes 

    real(dl), dimension(1:size(fld,dim=1),1:2) :: df_rel, df_tot
    real(dl) :: spec(1:size(fld,dim=1)/2+1)
    real(dl) :: norm, dk, keff
    integer :: i,j

    ! Need to remove the len ugliness
    norm = 1._dl / sqrt(2._dl*len*rho)
    dk = twopi/len

    df_rel = 0._dl
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl/sqrt(keff)
    enddo
    spec = spec * norm
    do i = 1,nFld; do j=1,2
       if (modes(2)) call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo; enddo

    df_tot = 0.
    spec = 0._dl
    do i=1,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl / (keff**2+m2eff)**0.25 
    enddo
    spec = spec * norm
    do i=1,nFld; do j=1,2
       if (modes(1)) call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo; enddo

    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo
  end subroutine initialise_fluctuations_kg_long

  
  subroutine initialise_fluctuations_bogoliubov_long(fld, rho, nCut)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), intent(in) :: rho
    integer, intent(in), optional :: nCut

    real(dl) :: df_rel(1:nLat,1:2), df_tot(1:nLat,1:2), spec(1:nLat/2+1)
    integer :: i,j
    real(dl) :: norm, dk, keff
    real(dl) :: lameff, nu_

    nu_ = nu  ! Fix this nonlocality
    lameff = del*sqrt(2._dl/nu_)
    norm = 1._dl / sqrt(2._dl*len*rho)
    dk = twopi / len   ! Fix this nonlocality
    
    ! Relative modes
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = sqrt( (keff**2+2._dl)/(keff*sqrt(keff**2+4._dl)) )
    enddo
    spec = spec * norm
    do j=1,2
       call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo
    
    ! Now get total modes
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = sqrt( (keff**2 + 2._dl - 4._dl*nu_)  / sqrt(keff**2+4._dl*nu_*(lameff**2-1._dl)) / sqrt(keff**2+4._dl-4._dl*nu_*(lameff**2+1._dl)) )
    enddo
    spec = spec * norm
    do j=1,2
       call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo
    
    ! Now remix the fields and add them to the background
    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo   
  end subroutine initialise_fluctuations_bogoliubov_long
  
end Module Fluctuations
