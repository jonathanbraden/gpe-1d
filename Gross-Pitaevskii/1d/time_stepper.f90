module TimeStepping
  use constants, only : dl, twopi
  use integrator
  use output

  implicit none
  
  type TimeStepper
     real(dl) :: tcur = 0.
     real(dl) :: dt = 1.
     real(dl) :: dt_out = 1.
     integer :: out_size = 1
     integer :: n_out_steps = 1
     integer :: n_steps = 1
  end type TimeStepper
  
contains

  subroutine time_evolve_stepper(fld, stepper, output_, log_)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(TimeStepper), intent(in) :: stepper
    logical, intent(in), optional :: output_, log_

    logical :: output, log
    integer :: i,j

    output = .true.; if (present(output_)) output=output_
    log = .true.; if (present(log_)) log=log_
    
    if (output) call output_fields_binary(fld)
    if (log) call output_log_file(fld, 0._dl)
    
    do i=1,stepper%n_out_steps
       do j=1,stepper%out_size
          call gl10(yvec, stepper%dt)
       enddo
       if (output) call output_fields_binary(fld)
       if (log) call output_log_file(fld, stepper%dt * i *stepper%out_size)
    enddo
  end subroutine time_evolve_stepper

  subroutine set_time_steps_oscillator(stepper, omega, w_samp)
    type(TimeStepper), intent(inout) :: stepper
    real(dl), intent(in) :: omega
    integer, intent(in) :: w_samp
    
    stepper%dt = (twopi/omega) / dble(w_samp)
    stepper%tcur = 0._dl

    ! out_size, n_out_steps

    stepper%dt_out = stepper%dt * stepper%out_size
    stepper%n_steps = stepper%out_size * stepper%n_out_steps
  end subroutine set_time_steps_oscillator

  subroutine set_time_steps_preheating(stepper)
    type(TimeStepper), intent(inout) :: stepper
  end subroutine set_time_steps_preheating
  
end module TimeStepping
