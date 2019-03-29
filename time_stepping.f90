module time_stepping

  use deriven, only: derivs
  use ghosts, only: set_ghosts

  implicit none

contains
  subroutine leap_frog_integrator(n_real, n_ghosts, n_bound, t, dt, pos, vel, a, mass, rho, sml, gamma, p, cs, u, du)
    integer, intent(in) :: n_real, n_ghosts, n_bound
    real, intent(inout) :: mass(n_real+n_ghosts), pos(n_real+n_ghosts), vel(n_real+n_ghosts), a(n_real+n_ghosts), &
    rho(n_real+n_ghosts), sml(n_real+n_ghosts), t,  u(n_real+n_ghosts), gamma, du(n_real+n_ghosts)
    real, intent(in) :: dt
    real, intent(out) :: p(n_real+n_ghosts), cs(n_real+n_ghosts)
    real :: a_0(n_real + n_ghosts), du_0(n_real + n_ghosts)

    !update what t we are on
    t = t+dt

    ! set arrays of the accel and du values used in the scheme final step.
    a_0(:) = a(:)
    du_0(:) = du(:)

    !probably call set_ghosts too many times.

    !Leap frog scheme according to q8 for velocity and internal energy, no pos equivalent so don't need to update.
    pos(:) = pos(:) + dt*vel(:) + 0.5*(dt)**(2)*a(:)
    call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)

    vel(:) = vel(:) + dt*a(:)
    u(:) = u(:) + dt*du(:)

    call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)


    call derivs(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, gamma, p, cs, u, du)

    vel(:) = vel(:) + 0.5*dt*(a(:) - a_0(:))
    u(:) = u(:) + 0.5*dt*(du(:) - du_0(:))
    call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)


  end subroutine leap_frog_integrator

  subroutine step_size(n_real, time_cond, sml, cs, dt)
    integer, intent(in) :: n_real
    real, intent(in) :: time_cond, sml(n_real), cs(n_real)
    real, intent(out) :: dt
    real :: time_array(n_real)

    !Find the minimum step size according to the courant condition
    !dt < h/cs, h and cs can vary so find minimum for all particles.
    time_array(:) = time_cond * sml(:)/cs(:)
    dt = minval(time_array)

  end subroutine step_size

end module time_stepping
