module bubble_extraction
  use constants, only : dl, twopi
  use fftw3
  implicit none
  
contains
  
  ! To do, adjust the lower bound here to depend on barrier location
  !
  ! Should tune thresholds to make use of some prescribed bubble radius, and field amplitude from instanton.
  !
  ! Known Bugs: Have to deal with periodic boundary, or else a bubble hitting boundary is counted twice 
 function count_bubbles(fld) result(nBub)
    real(dl), dimension(:), intent(in) :: fld
    integer, dimension(1:nt) :: nBub
    integer :: bSize  ! make this in an input parameter

    real(dl), dimension(1:size(fld)) :: cphi
    integer :: n, bSize
    integer, parameter :: nt = 4
    integer, dimension(1:nt) :: nv, nbound
    real(dl), dimension(1:nt) :: thresh
    logical, dimension(1:nt) :: boundary
    integer :: i,j
    
    n = size(fld); cphi = cos(fld)
    thresh = (/ 0., 0.25, 0.5, 0.75 /)
    nv = 0; nBub = 0; bSize = 15
    boundary = .false.; nbound = 0
    
    do j=1,nt
       if ( thresh(j) < cphi(1) ) then
          nv(j) = nv(j) + 1
          nbound(j) = nbound(j)+1
          boundary(j) = .true.
       endif
    enddo
    
    do i=2,n-1
       do j=1,nt
          if ( thresh(j) < cphi(i) ) then
             nv(j) = nv(j)+1
             if (boundary(j)) nbound(j) = nbound(j) + 1
          else
             if (nv(j) > bSize) nBub(j) = nBub(j) + 1
             nv(j) = 0
             if (boundary(j)) boundary(j) = .false.
          endif
       enddo
    enddo

    do j=1,nt
       if ( thresh(j) < cphi(n) ) then
          nv(j) = nv(j) + 1
          if (nbound(j) > 0) then
             if (nbound(j) <= bSize .and. nv(j)+nbound(j) > bSize) then
                nBub(j) = nBub(j) + 1
             endif
          endif
       endif
    enddo
    
  end function count_bubbles

  function int_stop(rFrac,dx,meff) result(i)
    real(dl), intent(in) :: rFrac, dx, meff
    integer :: i

    i = int(rFrac/(meff*dx))+1
  end function int_stop

  function count_bubbles_convolve(fld,scl) result(nBub)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), intent(in) :: scl
    integer :: nBub
    real(dl), dimension(1:size(fld)) :: f_sm

    ! Write this
  end function count_bubbles_convolve
  
  function count_bubbles_(fld,bSize) result(nBub)
    real(dl), dimension(:), intent(in) :: fld
    integer, intent(in) :: bSize
    
    real(dl), dimension(1:size(fld)) :: cphi
    integer :: n
    integer, parameter :: nt = 4
    integer, dimension(1:nt) :: nv, nbound, nBub
    real(dl), dimension(1:nt) :: thresh
    logical, dimension(1:nt) :: boundary
    integer :: i,j
    
    n = size(fld); cphi = cos(fld)
    thresh = (/ 0., 0.25, 0.5, 0.75 /)
    
    nv = 0; nBub = 0
    boundary = .false.; nbound = 0
    
    do j=1,nt
       if ( thresh(j) < cphi(1) ) then
          nv(j) = nv(j) + 1
          nbound(j) = nbound(j)+1
          boundary(j) = .true.
       endif
    enddo
    
    do i=2,n-1
       do j=1,nt
          if ( thresh(j) < cphi(i) ) then
             nv(j) = nv(j)+1
             if (boundary(j)) nbound(j) = nbound(j) + 1
          else
             if (nv(j) > bSize) nBub(j) = nBub(j) + 1
             nv(j) = 0
             !if (boundary(j)) boundary(j) = .false.  ! redundant
          endif
       enddo
    enddo

    do j=1,nt
       if ( thresh(j) < cphi(n) ) then
          nv(j) = nv(j) + 1
          if (nbound(j) > 0) then
             if (nbound(j) <= bSize .and. nv(j)+nbound(j) > bSize) then
                nBub(j) = nBub(j) + 1
             endif
          endif
       endif
    enddo
    
  end function count_bubbles_
 
  real(dl) function mean_cos(fld) result(cphi)
    real(dl), dimension(:), intent(in) :: fld
    
    cphi = sum(cos(fld))/dble(size(fld))
  end function mean_cos

  ! Add a smoothing here, then cosine, then clip
  function mean_clip(fld,th) result(mean)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), intent(in) :: th
    real(dl), dimension(1:size(fld)) :: cphi, clip
    real(dl) :: mean
    
    cphi = cos(fld)
    clip = 0._dl
    where (cphi > th) clip = 1._dl
    mean = sum(clip) / dble(size(fld))
  end function mean_clip

! Add in the smoothing window
  subroutine gaussian_smooth(fld_s,fld,tPair,kc)
    real(dl), dimension(:), intent(out) :: fld_s
    real(dl), dimension(:), intent(in) :: fld
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: kc
    integer :: i, nn
    
    nn = size(fld)/2+1
    tPair%realSpace(:) = fld(:)
    call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
    do i=1,nn
       tPair%specSpace(i) = tPair%specSpace(i)!*exp(-(i-1)**2*(kk/kc)**2)  ! Fill in filter
    enddo
    call fftw_execute_dft_c2r(tPair%planb,tPair%specSpace, tPair%realSpace)
    fld_s(:) = tPair%realSpace(:)
  end subroutine gaussian_smooth

  subroutine smooth_filter(fld_s,fld,wk,tPair)
    real(dl), dimension(:), intent(out) :: fld_s
    real(dl), dimension(:), intent(in) :: fld
    complex(dl), dimension(:), intent(in) :: wk
    type(transformPair1D), intent(inout) :: tPair

    tPair%realSpace(:) = fld(:)
    call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)
    tPair%specSpace(:) = tPair%specSpace(:)*conjg(wk(:))
    call fftw_execute_dft_c2r(tPair%planb,tPair%specSpace,tPair%realSpace)
    fld_s(:) = tPair%realSpace(:)
  end subroutine smooth_filter

  ! Fix this so it's the symmetric one
  subroutine tukey_filter(n,wk,alph,tPair)
    integer, intent(in) :: n
    complex(dl), dimension(1:n/2+1), intent(out) :: wk
    real(dl), intent(in) :: alph
    type(transformPair1D), intent(inout) :: tPair

    real(dl), dimension(1:n) :: f
    integer :: i, n_b, n_int
    real(dl) :: nt, pi
    pi = 0.5_dl*twopi
    nt = 0.5_dl*alph*n
    n_b = int(nt)

    do i=1,n_b-1
       f(i) = 0.5_dl*( 1._dl+cos(pi*(i/nt - 1._dl)) )
    enddo
    f(n_b:i-n_b) = 1._dl
    do i=i-n_b,n
       f(i) = 0.5_dl*( 1._dl*cos(pi*(i/nt - 2._dl/alph + 1._dl)) )
    enddo
  end subroutine tukey_filter
    
  !>@brief
  !> Apply a real-space top hat filter to the data
  !
  !>@todo Write the functionality, including periodicity
  subroutine box_filter_k(n,ns,wk,tPair)
    integer, intent(in) :: n, ns
    complex(dl), dimension(1:n/2+1), intent(out) :: wk
    type(transformPair1D), intent(inout) :: tPair

  end subroutine box_filter_k

  ! Produce a box filter smoothing over 2*ns+1 lattice sites
  subroutine box_filter_r(n,ns,w)
    integer, intent(in) :: n, ns
    real(dl), dimension(1:n), intent(out) :: w
    integer :: i; real(dl) :: norm
    ! Add some checks that ns is bigger than 0
    norm = 1./dble(2*ns+1)
    w = 0._dl
    w(1:ns+1) = norm; w(n-ns+1:n) = norm
  end subroutine box_filter_r
  
  subroutine make_fourier_window(w,wk,tPair)
    real(dl), dimension(:), intent(in) :: w
    complex(dl), dimension(:), intent(out) :: wk
    type(transformPair1D), intent(inout) :: tPair

    ! Normalize window
    tPair%realSpace(:) = w(:) / sum(w)
    call fftw_execute_dft_r2c(tPair%planf,tPair%realSpace,tPair%specSpace)    
    wk(:) = tPair%specSpace / size(w)
  end subroutine make_fourier_window
    
end module bubble_extraction
