module Output
  use constants, only : dl
  use utils, only : newunit
  use eom ! needed for model parameters, etc.  Refactor to remove
  implicit none

contains
  
  ! Debug this so that all the parameters are passed in as arguments
  subroutine write_header(oFile, dt)
    integer, intent(in) :: oFile
    real(dl), intent(in) :: dt

    write(oFile,*) "# GPE Simulation Lattice Parameters :"
    write(oFile,*) "# n_lattice = ",nLat," dx = ",dx
    write(oFile,*) "# Time Stepping Parameters :"
    write(oFile,*) "# dt = ",dt
    write(oFile,*) "# Model Parameters :"
    write(oFile,*) "# nu = ",nu,", g = ",gs," w = ",omega
    !write(oFile,*) "# rho = ",rho
    write(oFile,*) "# delta = ",del
    write(oFile,*) "#============================"
    write(oFile,*) "# t_{heal}   2\sqrt{nu}t_{heal}   rho"
  end subroutine write_header
  
  subroutine output_log_file(fld, t, dt, fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    real(dl), intent(in) :: t, dt
    character(80), intent(in), optional :: fName

    character(80) :: fn
    integer, save :: oFile
    logical :: o

    fn = 'log.out'
    if (present(fName)) fn = trim(fName)

    inquire(file=trim(fn), opened=o)
    if (.not.o) then
       open(unit=newunit(oFile), file=trim(fn))
       call write_header(oFile, dt)
    endif
    write(oFile,*) t, 2.*sqrt(nu)*t, sum(fld**2)/size(fld,dim=1)
    
  end subroutine output_log_file
  
  subroutine output_fields_binary(fld, fName)
    real(dl), dimension(:,:,:), intent(in) :: fld
    character(80), intent(in), optional :: fName

    logical :: o
    character(80) :: fn
    integer, save :: oFile

    fn = 'fields.bin'
    if (present(fName)) fn = trim(fName)

    inquire(file=trim(fn), opened=o)
    if (.not.o) open(unit=newunit(oFile), file=fn, access='stream')

    write(oFile) fld
  end subroutine output_fields_binary
    
end module Output
