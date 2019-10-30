module bubble_extraction

  use constants, only : dl
  implicit none
  
contains

  ! To do, adjust the lower bound here to depend on barrier location
  !
  ! Should tune thresholds to make use of some prescribed bubble radius, and field amplitude from instanton.
  !
  ! Known Bugs: Have to deal with periodic boundary, or else a bubble hitting boundary is counted twice
  function count_bubbles(fld) result(nBub)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), dimension(1:size(fld)) :: cphi
    integer :: n, bSize
    integer, parameter :: nt = 4
    integer, dimension(1:nt) :: nv, nbound, nBub
    real(dl), dimension(1:nt) :: thresh
    logical, dimension(1:nt) :: boundary
    integer :: i,j
    
    n = size(fld); cphi = cos(fld)
    thresh = (/ 0., 0.25, 0.5, 0.75 /)
    nv = 0; nBub = 0; bSize = 5
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

  real(dl) function mean_cos(fld) result(cphi)
    real(dl), dimension(:), intent(in) :: fld
    real(dl), dimension(1:size(fld)) :: cp
    
    cphi = sum(cos(fld))/dble(size(fld))
  end function mean_cos
    
  subroutine minkowski_functional(thresh)
    real(dl), dimension(:), intent(in) :: thresh
  end subroutine minkowski_functional
  
end module bubble_extraction
