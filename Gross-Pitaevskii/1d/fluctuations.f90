! TO DO : Remove the len nonlocality (and any others) in this code
! TO DO : Finish the wrapper function

module Fluctuations
  use constants, only : dl, twopi
  use gaussianRandomField
  use gpe_model ! remove this dependence. It stores grid params, etc.
  
  implicit none

  type SpecParams
     real(dl) :: cos_phi, dk
     real(dl) :: m2eff, nu, lamEff
     integer :: num_atoms
     integer :: nCut = 1
     character(8) :: type = 'BOGO'
     logical, dimension(2) :: modes = (/.true.,.true./) ! Update for nFlds
     real(dl) :: rho, len
  end type SpecParams
  
contains

  ! Get this to take in model parameters
  function make_spec_params(rho, len, nu, lamEff, type, modes, nCut, fv_) result(params)
    real(dl), intent(in) :: rho, len, nu, lamEff
    character(*), intent(in) :: type
    logical, dimension(1:2), intent(in) :: modes
    integer, intent(in) :: nCut
    logical, intent(in), optional :: fv_  ! Don't make optional
    
    type(SpecParams) :: params

    logical :: fv
    integer :: sig
    
    fv = .true.; if (present(fv_)) fv = fv_
    params%cos_phi = 1._dl; if (fv) params%cos_phi = -1._dl

    params%nCut = nCut
    params%rho = rho  ! Merge w/ len into num_atoms
    params%len = len  ! Merge w/ len into num_atoms
    params%nu = nu
    params%lamEff = lamEff  ! Add appropriate conversions
    params%m2eff = 4._dl*nu*(lamEff**2 + params%cos_phi*1._dl) 
    params%type = trim(type)
    params%modes(1:2) = modes(1:2)
    params%dk = twopi / len
  end function make_spec_params

  subroutine print_spec_params(params)
    type(SpecParams), intent(in) :: params

    print*,"=========================================="
    print*,"Fluctuation Properties"
    print*,"------------------------------------------"
    print*,"type of fluctuations   : ", params%type
    print*,"m^2_eff                : ", params%m2eff
    print*,"Initialise total modes : ", params%modes(1)
    print*,"Initialise rel modes   : ", params%modes(2)
    print*,"Spectral cutff number  : ", params%nCut
    print*,"Number of atoms        : ", params%len * params%rho
    print*,"Cos(\phi_bg)           : ", params%cos_phi
    print*,"==========================================="
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
    case ('KG_R')
       call initialise_fluctuations_kg_real(fld, params)
    case default
       call initialise_fluctuations_white(fld, params)
    end select
  end subroutine initialise_fluctuations

  subroutine initialise_fluctuations_bogoliubov(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    complex(dl), dimension(1:size(fld,dim=1)) :: df_pos, df_neg
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_neg, spec_pos  ! Adjust size
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_u, spec_v      ! Adjust size
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: keff               ! Adjust size
    integer :: i,j
    real(dl) :: norm, dk
    real(dl) :: lameff, nu_
    real(dl) :: m2eff
    real(dl) :: num_atoms
    integer :: nCut

    nCut = params%nCut
    lameff = params%lamEff
    nu_ = params%cos_phi * params%nu  ! Make sure I sync this sign correctly below
    m2eff = params%m2eff
    num_atoms = params%num_atoms
    
    norm = 1._dl / sqrt(2._dl*num_atoms)  ! check factor of 2
    dk = params%dk
    keff = (/ (dk*(i-1), i=1,size(spec_pos)) /)

    spec_pos = 0._dl; spec_neg = 0._dl
    if (params%cos_phi < 0.) then
       spec_neg(2:) = sqrt( (keff(2:)**2+2._dl) / (keff(2:)*sqrt(keff(2:)**2+4._dl)) )
       spec_pos(2:) = sqrt( (keff(2:)**2 + 2._dl + 4._dl*nu_)  / sqrt(keff(2:)**2+m2eff) / sqrt(keff(2:)**2+4._dl+4._dl*nu_*(lameff**2+1._dl)) )
    else
       spec_neg(2:) = sqrt( (keff(2:)**2 + 2._dl + 4._dl*nu_)  / sqrt(keff(2:)**2+m2eff) / sqrt(keff(2:)**2+4._dl+4._dl*nu_*(lameff**2+1._dl)) )
       spec_pos(2:) = sqrt( (keff(2:)**2+2._dl) / (keff(2:)*sqrt(keff(2:)**2+4._dl)) )
    endif

    df_pos = 0._dl
    if (params%modes(1)) then
       call convert_spec_psi_to_u_and_v(spec_pos, spec_u, spec_v, noise_floor_=1._dl)
       call generate_1dGRF_cmplx(df_pos, spec_u, spec_v)
       df_pos = df_pos * norm
    endif
    
    df_neg = 0._dl
    if (params%modes(2)) then
       call convert_spec_psi_to_u_and_v(spec_neg, spec_u, spec_v, noise_floor_=1._dl)
       call generate_1dGRF_cmplx(df_neg, spec_u, spec_v)
       df_neg = df_neg * norm
    endif
       
    fld(:,1,1) = fld(:,1,1) + sqrt(0.5_dl)*( real(df_pos) + real(df_neg) )
    fld(:,1,2) = fld(:,1,2) + sqrt(0.5_dl)*( real(df_pos) - real(df_neg) )
    fld(:,2,1) = fld(:,2,1) + sqrt(0.5_dl)*( aimag(df_pos) + aimag(df_neg) )
    fld(:,2,2) = fld(:,2,2) + sqrt(0.5_dl)*( aimag(df_pos) - aimag(df_neg) )
   
  end subroutine initialise_fluctuations_bogoliubov

  subroutine initialise_fluctuations_bogoliubov_real_imag(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    ! Write this to initialise real and imaginary modes directly
  end subroutine initialise_fluctuations_bogoliubov_real_imag
  
  subroutine initialise_fluctuations_white(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    complex(dl), dimension( 1:size(fld,dim=1) ) :: df
    real(dl), dimension( 1:size(fld,dim=1)/2+1 ) :: spec
    integer :: num_fld, n
    integer :: i, j

    real(dl) :: norm, num_part
    integer :: nCut

    nCut = params%nCut
    num_part = params%rho * params%len
    
    norm = 1._dl / sqrt(2._dl*num_part)  ! Check factors of 2 in here
    
    n = size(fld,dim=1)
    num_fld = size(fld,dim=3)
    
    spec = 0._dl
    do i=2,n/2; do j=1,2
       spec(i) = 1._dl
    enddo; enddo

    do i = 1,num_fld
       !call generate_1dGRF_cmplx(df, spec(1:nCut))
       df = df * norm
       fld(:,1,i) = fld(:,1,i) + real(df)  ! Check unmixing factor
       fld(:,2,i) = fld(:,2,i) + aimag(df) ! Check unmixing factor
    enddo
  end subroutine initialise_fluctuations_white

  ! This if broken.  Figure out how to fix
  ! Fix the nCut stuff (check with random_fields file)
  ! Fix normalisation conventions between different subroutins
  ! Modularise the mixing back into psi_1 and psi_2
  subroutine initialise_fluctuations_kg(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams) :: params

    complex(dl), dimension(1:size(fld,dim=1)) :: df_pos, df_neg
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_pos, spec_neg, keff
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_u, spec_v
    real(dl) :: m2eff, dk, norm
    integer :: nCut, i,j
    real(dl) :: num_atoms

    m2eff = params%m2eff
    nCut = params%nCut
    
    norm = 1._dl / sqrt(2._dl*params%num_atoms) 
    dk = params%dk
    keff = (/ ((i-1)*dk, i=1,size(keff))  /)

    spec_pos = 0.; spec_neg = 0.
    if (params%cos_phi < 0.) then
       spec_pos(2:) = 1._dl/(keff(2:)**2+m2eff)**0.25
       spec_neg(2:) = 1._dl/sqrt(keff(2:))  
    else
       spec_pos(2:) = 1._dl/sqrt(keff(2:))
       spec_neg(2:) = 1._dl/(keff(2:)**2+m2eff)**0.25
    endif

    df_pos = 0._dl
    if (params%modes(1)) then
       call convert_spec_psi_to_u_and_v(spec_pos, spec_u, spec_v, noise_floor_=0.)
       call generate_1dGRF_cmplx( df_pos, spec_u, spec_v )
       df_pos = df_pos*norm
    endif

    df_neg = 0._dl
    if (params%modes(2)) then
       call convert_spec_psi_to_u_and_v(spec_neg, spec_u, spec_v, noise_floor_=0.)
       call generate_1dGRF_cmplx( df_neg, spec_u, spec_v )
       df_neg = df_neg*norm
    endif
    
    fld(:,1,1) = fld(:,1,1) + sqrt(0.5_dl)*( real(df_pos) + real(df_neg) ) 
    fld(:,2,1) = fld(:,2,1) + sqrt(0.5_dl)*( aimag(df_pos) - aimag(df_neg) )
    fld(:,1,1) = fld(:,1,1) + sqrt(0.5_dl)*( real(df_pos) + real(df_neg) ) 
    fld(:,2,2) = fld(:,2,2) + sqrt(0.5_dl)*( aimag(df_pos) - aimag(df_neg) )
  end subroutine initialise_fluctuations_kg
  
  subroutine initialise_fluctuations_kg_real(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams) :: params

    real(dl), dimension(1:size(fld,dim=1),1:2) :: df_rel, df_tot
    real(dl) :: spec(1:size(fld,dim=1)/2+1)
    real(dl) :: m2eff, dk, keff, norm
    integer :: nCut, i,j
    logical, dimension(2) :: modes
    real(dl) :: len, num_atoms

    len = params%len
    num_atoms = params%len * params%rho
    
    modes(1:2) = params%modes(1:2)    
    m2eff = params%m2eff
    nCut = params%nCut
    
    norm = 1._dl / sqrt(2._dl*num_atoms)
    dk = twopi/len

    df_rel = 0._dl
    spec = 0._dl
    do i=2,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl/sqrt(keff)
    enddo
    spec = spec * norm
    do j=1,2
       if (modes(2)) call generate_1dGRF(df_rel(:,j), spec(1:nCut), .false.)
    enddo

    df_tot = 0.
    spec = 0._dl
    do i=1,nLat/2
       keff = (i-1)*dk
       spec(i) = 1._dl / (keff**2+m2eff)**0.25 
    enddo
    spec = spec * norm
    do j=1,2
       if (modes(1)) call generate_1dGRF(df_tot(:,j), spec(1:nCut), .false.)
    enddo

    do j=1,2
       fld(:,j,1) = fld(:,j,1) + sqrt(0.5_dl)*( df_tot(:,j) + df_rel(:,j) )
       fld(:,j,2) = fld(:,j,2) + sqrt(0.5_dl)*( df_tot(:,j) - df_rel(:,j) )
    enddo
  end subroutine initialise_fluctuations_kg_real
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More verbose versions of calls, with parameters explicitly listed
! This will be deleted eventually.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  subroutine initialise_relative_phase_fluctuations(fld, spec, f1, f2)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), dimension(:), intent(in) :: spec
    integer, intent(in) :: f1, f2

    real(dl), dimension(1:size(fld,dim=1)) :: dphase
    real(dl), dimension(1:size(fld,dim=1)) :: rho

    dphase = 0._dl  ! Delete this once I call the GRF
    ! Insert code to generate phase fluctuations using GRF generator

    ! Extract current values of phase
  end subroutine initialise_relative_phase_fluctuations

  subroutine initialise_relative_density_fluctuations(fld, spec, f1, f2, phi0)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), dimension(:), intent(in) :: spec
    integer, intent(in) :: f1, f2
    real(dl), intent(in) :: phi0
    
    real(dl), dimension(1:size(fld,dim=1)) :: drho
    real(dl) :: rho_cur

    rho_cur = 1._dl
    drho = 0._dl
    ! Insert code to generate density fluctuations
    fld(:,1,f2) = sqrt(rho_cur+0.5_dl*drho)*cos(0.5_dl*phi0)
    fld(:,2,f2) = sqrt(rho_cur+0.5_dl*drho)*sin(0.5_dl*phi0)
    fld(:,1,f1) = sqrt(rho_cur-0.5_dl*drho)*cos(0.5_dl*phi0)
    fld(:,2,f1) = sqrt(rho_cur-0.5_dl*drho)*sin(0.5_dl*phi0)
  end subroutine initialise_relative_density_fluctuations

  !>@brief
  !> Converts a spectrum for $\psi^2$ into separate u and v spectra.
  !> Assumes the "quantum noise" lower bound on the spectrum is 1/2
  subroutine convert_spec_psi_to_u_and_v(spec, spec_u, spec_v, noise_floor_)
    real(dl), dimension(:), intent(in) :: spec
    real(dl), dimension(1:size(spec)), intent(out) :: spec_u, spec_v
    real(dl), intent(in), optional :: noise_floor_

    real(dl) :: noise_floor

    noise_floor = 0.5_dl; if (present(noise_floor_)) noise_floor=noise_floor_
    
    spec_u = sqrt(spec**2 + noise_floor)  ! Add some checks that this is positive
    spec_v = sqrt(spec**2 - noise_floor)  ! Add some checks that this is positive
  end subroutine convert_spec_psi_to_u_and_v

  ! Debug this to make sure the overall normalization is correct
  ! and that I've used the correct signs
  subroutine convert_spec_psi_to_real_and_imag(spec, spec_r, spec_i, noise_floor_)
    real(dl), dimension(:), intent(in) :: spec
    real(dl), dimension(1:size(spec)), intent(out) :: spec_r, spec_i
    real(dl), intent(in), optional :: noise_floor_

    real(dl) :: noise_floor

    noise_floor = 0.5_dl; if (present(noise_floor_)) noise_floor=noise_floor_

    ! Check normalisation on these, and also sign
    spec_r = sqrt(spec**2 + noise_floor) - sqrt(spec**2-noise_floor)
    spec_i = sqrt(spec**2 + noise_floor) + sqrt(spec**2-noise_floor) 
  end subroutine convert_spec_psi_to_real_and_imag
  
end Module Fluctuations
