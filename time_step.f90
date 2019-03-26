module time_step

  implicit none

contains
  subroutine step_size(n_real, time_cond, sml, cs, dt)
    integer, intent(in) :: n_real
    real, intent(in) :: time_cond, sml(n_real), cs(n_real)
    real, intent(out) :: dt
    real :: time_array(n_real)

    time_array(:) = time_cond * sml(:)/cs(:)
    dt = minval(time_array)

  end subroutine step_size

end module time_step
