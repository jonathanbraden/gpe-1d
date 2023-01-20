module TimeStepping
  use constants, only : dl, twopi
  use integrator
  use output

  implicit none
  
  type TimeStepper
     real(dl) :: tcur = 0.
     real(dl) :: dt = 1.
     integer :: out_size = 1
     integer :: n_out_steps = 1
  end type TimeStepper
  
contains

  subroutine time_evolve_stepper(fld, stepper, output_, log_, verbose_)
    real(dl), dimension(:,:,:), intent(inout) :: fld
    type(TimeStepper), intent(in) :: stepper
    logical, intent(in), optional :: output_, log_, verbose_

    logical :: output, log, verbose
    integer :: i,j
    real(dl) :: dt_out

    dt_out = stepper%dt * stepper%out_size

    output = .true.; if (present(output_)) output=output_
    log = .true.; if (present(log_)) log=log_
    verbose = .false.; if (present(verbose_)) verbose = verbose_
    
    if (output) call output_fields_binary(fld)
    if (log) call output_log_file(fld, 0._dl)
    
    do i=1,stepper%n_out_steps
       if (verbose) print*,"Starting output step ",i
       do j=1,stepper%out_size
          call gl10(yvec, stepper%dt)
       enddo
       if (output) call output_fields_binary(fld)
       if (log) call output_log_file(fld, i*dt_out)
    enddo
  end subroutine time_evolve_stepper

  subroutine set_time_steps_oscillator(stepper, omega, w_samp, out_size, t_final)
    type(TimeStepper), intent(inout) :: stepper
    real(dl), intent(in) :: omega
    integer, intent(in) :: w_samp, out_size
    real(dl), intent(in) :: t_final
    real(dl) :: dt_out

    stepper%tcur = 0._dl
    stepper%dt = (twopi/omega) / dble(w_samp)
    stepper%out_size = out_size

    dt_out = stepper%dt * stepper%out_size
    stepper%n_out_steps = int(t_final / dt_out)

  end subroutine set_time_steps_oscillator

  subroutine set_time_steps_preheating(stepper)
    type(TimeStepper), intent(inout) :: stepper
  end subroutine set_time_steps_preheating

  subroutine print_time_stepper(stepper)
    type(TimeStepper), intent(in) :: stepper

    print*,"========================"
    print*,"TimeStepper properties :"
    print*,"------------------------"
    print*,"dt : ",stepper%dt
    print*,"dt_out : ",stepper%dt * stepper%out_size
    print*,"time steps per output step (dt_out/dt) : ",stepper%out_size
    print*,"total number of output steps : ",stepper%n_out_steps
    print*,"total number of time steps : ",stepper%out_size * stepper%n_out_steps
    print*,"final output time : ",stepper%out_size * stepper%n_out_steps * stepper%dt
    print*,"tcur : ",stepper%tcur
    print*,"========================"
  end subroutine print_time_stepper
  
end module TimeStepping
