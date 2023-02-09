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
     real(dl) :: num_atoms
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
    params%num_atoms = rho*len
  end function make_spec_params

  subroutine print_spec_params(params)
    type(SpecParams), intent(in) :: params

    print*,""
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
    print*,""
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
    case ('BOGO_DP')
       call initialise_fluctuations_phase_and_density_bogoliubov(fld, params, 1._dl, 0._dl)
    case ('KG_DP')
       call initialise_fluctuations_phase_and_density_kg(fld, params, 1._dl, 0._dl)
    case default
       call initialise_fluctuations_bogoliubov(fld, params)
    end select
  end subroutine initialise_fluctuations

  ! Debug this, then replace code below
  subroutine spectrum_bogoliubov(spec, params, is_relative)
    real(dl), dimension(:), intent(out) :: spec
    type(SpecParams), intent(in) :: params
    logical, intent(in) :: is_relative

    real(dl), dimension(1:size(spec)) :: keff
    real(dl) :: dk, norm
    real(dl) :: nu_, m2eff, lameff
    integer :: i, nCut

    nu_ = params%cos_phi * params%nu
    lameff = params%lamEff
    m2eff = params%m2eff
    
    dk = params%dk
    norm = 1./sqrt(2.*params%num_atoms)
    nCut = params%nCut
    ! Add a check on size of nCut
    
    keff = (/ ((i-1)*dk, i=1,size(keff)) /)
    spec = 0._dl
    if (is_relative) then
       spec(2:nCut) = sqrt( (keff(2:nCut)**2+2._dl+4._dl*nu_) / sqrt(keff(2:nCut)**2+m2eff) / sqrt(keff(2:nCut)**2+4._dl+4._dl*nu_*(lameff**2+1._dl)) )  ! Check this one is correct
    else
       spec(2:nCut) = sqrt( (keff(2:nCut)**2+2._dl) / ( keff(2:nCut)*sqrt(keff(2:nCut)**2+4._dl) ) )
    endif
  end subroutine spectrum_bogoliubov

  subroutine spectrum_kg(spec, params, is_relative)
    real(dl), dimension(:), intent(out) :: spec
    type(SpecParams), intent(in) :: params
    logical, intent(in) :: is_relative

    real(dl), dimension(1:size(spec)) :: keff
    real(dl) :: dk, norm
    real(dl) :: m2eff
    integer :: i, nCut

    nCut = params%nCut
    dk = params%dk
    norm = 1._dl / sqrt(2.*params%num_atoms)
    
    keff = (/ ((i-1)*dk, i=1,size(keff)) /)
    spec = 0._dl
    if (is_relative) then
       spec(2:nCut) = 1._dl / sqrt( keff(2:nCut)**2+m2eff)
    else
       spec(2:nCut) = 1._dl / sqrt(keff(2:nCut))
    endif
  end subroutine spectrum_kg

  ! Add subroutine to take the k->0 limit, including the overall normalization different
  
  subroutine initialise_fluctuations_bogoliubov(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params

    complex(dl), dimension(1:size(fld,dim=1)) :: df_pos, df_neg
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_neg, spec_pos  ! Adjust size
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_tot, spec_rel
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: spec_u, spec_v      ! Adjust size
    real(dl), dimension(1:size(fld,dim=1)/2+1) :: keff               ! Adjust size
    integer :: i,j
    real(dl) :: norm, dk
    real(dl) :: lameff, nu_
    real(dl) :: m2eff
    real(dl) :: num_atoms
    integer :: nCut

    nCut = params%nCut  ! Could incorporate this into the subroutine variable defs
    lameff = params%lamEff
    nu_ = params%cos_phi * params%nu  ! Make sure I sync this sign correctly below
    m2eff = params%m2eff
    num_atoms = params%num_atoms
    
    norm = 1._dl / sqrt(2._dl*num_atoms)  ! check factor of 2 (coincides w/ noise floor choice)
    dk = params%dk
    keff = (/ (dk*(i-1), i=1,size(spec_pos)) /)

    spec_tot = 0._dl; spec_rel = 0._dl
    spec_tot = sqrt( (keff(2:)**2+2._dl) / (keff(2:)*sqrt(keff(2:)**2+4._dl)) )
    spec_rel = sqrt( (keff(2:)**2+2._dl+4._dl*nu_) / sqrt(keff(2:)**2+m2eff) / sqrt(keff(2:)**2+4._dl+4._dl*nu_*(lameff**2+1._dl)) )

    if (params%cos_phi < 0.) then
       spec_neg(:) = spec_tot(:)
       spec_pos(:) = spec_rel(:)
    else
       spec_pos(:) = spec_tot(:) 
       spec_neg(:) = spec_rel(:)
    endif

    df_pos = 0._dl
    if (params%modes(1)) then
       call convert_spec_psi_to_u_and_v(spec_pos, spec_u, spec_v, noise_floor_=1._dl)
       call generate_1dGRF_cmplx(df_pos, spec_u(:nCut), spec_v(:nCut))
       df_pos = df_pos * norm
    endif

    df_neg = 0._dl
    if (params%modes(2)) then
       call convert_spec_psi_to_u_and_v(spec_neg, spec_u, spec_v, noise_floor_=1._dl)
       call generate_1dGRF_cmplx(df_neg, spec_u(:nCut), spec_v(:nCut))
       df_neg = df_neg * norm
    endif

    ! Patch these up to use cos_phi for the sign
    ! This will simplify the logic above
    fld(:,1,1) = fld(:,1,1) + sqrt(0.5_dl)*( real(df_pos) + real(df_neg) )
    fld(:,1,2) = fld(:,1,2) + sqrt(0.5_dl)*( real(df_pos) - real(df_neg) )
    fld(:,2,1) = fld(:,2,1) + sqrt(0.5_dl)*( aimag(df_pos) + aimag(df_neg) )
    fld(:,2,2) = fld(:,2,2) + sqrt(0.5_dl)*( aimag(df_pos) - aimag(df_neg) )
   
  end subroutine initialise_fluctuations_bogoliubov

  ! This one is redundant.  The KG one isn't
  ! However, I should write this regardless
  subroutine initialise_fluctuations_bogoliubov_real_imag(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params   
  end subroutine initialise_fluctuations_bogoliubov_real_imag

  ! Fix this now that it's broken
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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! More verbose versions of calls, with parameters explicitly listed
! This will be deleted eventually.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Debug this and clean it up
  subroutine initialise_fluctuations_phase_and_density(fld, spec_rel, spec_tot, norm, nf, phi0, rho0)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    real(dl), dimension(:), intent(in) :: spec_rel, spec_tot
    real(dl), intent(in) :: norm, nf, rho0, phi0

    real(dl), dimension(1:size(fld,1)) :: drho_tot, dphase_tot, drho_rel, dphase_rel
    real(dl), dimension(1:size(spec_rel)) :: spec_rho, spec_phase
    
    call convert_spec_psi_to_real_and_imag(spec_tot, spec_rho, spec_phase, noise_floor_=nf)
    call generate_1dGRF(drho_tot, spec_rho, .false.)
    call generate_1dGRF(dphase_tot, spec_phase, .false.)
    drho_tot = drho_tot*norm
    dphase_tot = dphase_tot*norm
    
    call convert_spec_psi_to_real_and_imag(spec_rel, spec_rho, spec_phase, noise_floor_=nf)
    call generate_1dGRF(drho_rel, spec_rho, .false.)
    call generate_1dGRF(dphase_rel, spec_phase, .false.)
    drho_rel = drho_rel*norm
    dphase_rel = dphase_rel*norm
    
    fld(:,1,1) = sqrt(1.+drho_tot-drho_rel)*cos(0.5*phi0 + 0.5_dl*(dphase_tot - dphase_rel) )
    fld(:,2,1) = sqrt(1.+drho_tot-drho_rel)*sin(-0.5*phi0 + 0.5_dl*(dphase_tot - dphase_rel) )
    fld(:,1,2) = sqrt(1.+drho_tot+drho_rel)*cos(0.5*phi0 + 0.5_dl*(dphase_tot + dphase_rel) )  ! check signs
    fld(:,2,2) = sqrt(1.+drho_tot+drho_rel)*sin(0.5*phi0 + 0.5_dl*(dphase_tot + dphase_rel) ) ! check signs
  end subroutine initialise_fluctuations_phase_and_density
  
  subroutine initialise_fluctuations_phase_and_density_bogoliubov(fld, params, rho0, phi0)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params
    real(dl), intent(in) :: rho0, phi0

    real(dl), dimension(1:size(fld,1)/2+1) :: spec_tot, spec_rel, spec_rho, spec_phase
    real(dl), dimension(1:size(fld,1)/2+1) :: keff
    real(dl), dimension(1:size(fld,1)) :: drho_tot, dphase_tot, drho_rel, dphase_rel
    real(dl) :: norm, dk
    real(dl) :: m2eff, lameff, nu_
    integer :: i


    ! I haven't put the norm in here correctly yet.
    dk = params%dk
    norm = 1._dl / sqrt(2._dl*params%num_atoms)
    keff = (/ ((i-1)*dk, i=1,size(keff)) /)

    nu_ = params%cos_phi * params%nu
    m2eff = params%m2eff
    lameff = params%lameff

    ! Add nCut, etc. in here.
    ! Even better, write external subroutine to calculate the spectrum, since I do it multiple times
    spec_rel = 0._dl
    spec_tot = 0._dl
    
    spec_tot(2:) = sqrt( (keff(2:)**2+2._dl) / ( keff(2:)*sqrt(keff(2:)**2+4._dl) ) )
    spec_rel(2:) = sqrt( (keff(2:)**2 + 2._dl + 4._dl*nu_) / sqrt(keff(2:)**2+m2eff) / sqrt(keff(2:)**2+4._dl+4._dl*nu_*(lameff**2+1._dl)) )

    if (params%modes(1)) then
       call convert_spec_psi_to_real_and_imag(spec_tot, spec_rho, spec_phase, noise_floor_=1.) ! Check noise floor norm
       call generate_1dGRF(drho_tot, spec_rho, .false.)
       call generate_1dGRF(dphase_tot, spec_phase, .false.)
    else
       drho_tot = 0._dl
       dphase_tot = 0._dl
    endif

    if (params%modes(2)) then
       call convert_spec_psi_to_real_and_imag(spec_rel, spec_rho, spec_phase, noise_floor_=1.) ! Check noise floor norm
       call generate_1dGRF(drho_rel, spec_rho, .false.)
       call generate_1dGRF(dphase_rel, spec_phase, .false.)
    else
       drho_rel = 0._dl
       dphase_rel = 0._dl
    endif
    
    drho_tot = drho_tot*norm
    dphase_tot = dphase_tot*norm
    drho_rel = drho_rel*norm
    dphase_rel = dphase_rel*norm

    ! An option here is to subtract off the mean field (if it exists), then add it back on later?
    fld(:,1,1) = sqrt(1.+drho_tot-drho_rel)*cos(0.5*phi0 + 0.5_dl*(dphase_tot - dphase_rel) ) 
    fld(:,2,1) = sqrt(1.+drho_tot-drho_rel)*sin(-0.5*phi0 + 0.5_dl*(dphase_tot - dphase_rel) )
    fld(:,1,2) = sqrt(1.+drho_tot+drho_rel)*cos(0.5*phi0 + 0.5_dl*(dphase_tot + dphase_rel) )  ! check signs
    fld(:,2,2) = sqrt(1.+drho_tot+drho_rel)*sin(0.5*phi0 + 0.5_dl*(dphase_tot + dphase_rel) ) ! check signs

  end subroutine initialise_fluctuations_phase_and_density_bogoliubov

  ! Copy the Bogoliubov above for the Klein-Gordon, but put in spectra explicitly so I can cut them off at high-k
  subroutine initialise_fluctuations_phase_and_density_kg(fld, params, rho0, phi0)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params
    real(dl), intent(in) :: rho0, phi0

    real(dl), dimension(1:size(fld,1)/2+1) :: spec_rho, spec_phase
    real(dl), dimension(1:size(fld,1)/2+1) :: keff
    real(dl), dimension(1:size(fld,1)) :: drho_tot, dphase_tot, drho_rel, dphase_rel
    real(dl) :: norm, dk
    real(dl) :: m2eff, nu_
    integer :: i

    dk = params%dk
    norm = 1._dl / sqrt(2._dl*params%num_atoms)
    keff = (/ ((i-1)*dk, i=1,size(keff)) /)

    nu_ = params%cos_phi * params%nu
    m2eff = params%m2eff

    if (params%modes(1)) then
       spec_rho = 0.
       spec_rho(2:) = keff(2:)
       spec_phase = 0.
       spec_phase(2:) = 1._dl/keff(2:)
       call generate_1dGRF(drho_tot, spec_rho, .false.)
       call generate_1dGRF(dphase_tot, spec_phase, .false.)
    else
       drho_tot = 0._dl
       dphase_tot = 0._dl
    endif

    ! Check these
    ! Also, compare normalization between different sims
    if (params%modes(2)) then
       spec_rho = 0._dl
       spec_rho(2:) = sqrt(keff(2:)**2+m2eff)
       spec_phase = 0._dl
       spec_phase(2:) = 1./sqrt(keff(2:)**2+m2eff)
       call generate_1dGRF(drho_rel, spec_rho, .false.)
       call generate_1dGRF(dphase_rel, spec_phase, .false.)
    else
       drho_rel = 0._dl
       dphase_rel = 0._dl
    endif

    drho_tot = drho_tot*norm
    dphase_tot = dphase_tot*norm
    drho_rel = drho_rel*norm
    dphase_rel = dphase_rel*norm

    ! An option here is to subtract off the mean field (if it exists), then add it back on later?
    fld(:,1,1) = sqrt(1.+drho_tot-drho_rel)*cos(0.5*phi0 + 0.5_dl*(dphase_tot - dphase_rel) ) 
    fld(:,2,1) = sqrt(1.+drho_tot-drho_rel)*sin(-0.5*phi0 + 0.5_dl*(dphase_tot - dphase_rel) )
    fld(:,1,2) = sqrt(1.+drho_tot+drho_rel)*cos(0.5*phi0 + 0.5_dl*(dphase_tot + dphase_rel) )  ! check signs
    fld(:,2,2) = sqrt(1.+drho_tot+drho_rel)*sin(0.5*phi0 + 0.5_dl*(dphase_tot + dphase_rel) ) ! check signs
    
  end subroutine initialise_fluctuations_phase_and_density_kg

  subroutine initialise_fluctuations_real_and_imaginary_bogoliubov(fld, params)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(SpecParams), intent(in) :: params
    
    real(dl), dimension(1:params%nCut) :: spec, spec_re, spec_im
    real(dl), dimension(1:size(fld,1)) :: psi_tot_re, psi_tot_im, psi_rel_re, psi_rel_im
    
    ! Use call to Bogoliubov spectrum here

    if ( params%modes(1) ) then
       ! Call the spectrum here
       spec_re = 0._dl
       spec_im = 0._dl
       call generate_1dGRF(psi_tot_re, spec_re, .false.)
       call generate_1dGRF(psi_tot_im, spec_im, .false.)
    endif
    
    if ( params%modes(2) ) then
       ! Call to get spectrum
       spec_re = 0._dl
       spec_im = 0._dl
       call generate_1dGRF(psi_rel_re, spec_re, .false.)
       call generate_1dGRF(psi_rel_im, spec_im, .false.)
    endif    
  end subroutine initialise_fluctuations_real_and_imaginary_bogoliubov
  
  
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

    !print*,"Noise floor is ",noise_floor
    ! Check normalisation on these, and also sign
    spec_r = 0._dl; spec_i = 0._dl

    ! Improve this to explicitly ignore the zero mode
    where (spec**2>noise_floor)
       spec_r = sqrt(spec**2 + noise_floor) - sqrt(spec**2-noise_floor)
       spec_i = sqrt(spec**2 + noise_floor) + sqrt(spec**2-noise_floor)
    end where
  end subroutine convert_spec_psi_to_real_and_imag
  
end Module Fluctuations
