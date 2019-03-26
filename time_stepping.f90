module time_stepping

  use deriven, only: derivs
  use ghosts, only: set_ghosts
  use output, only: write_output

  implicit none

contains
  subroutine leap_frog_integrator(n_real, n_ghosts, t, dx, dt, pos, vel, a, mass, rho, sml, gamma, p, cs, u, nfile)
    integer, intent(in) :: n_real, n_ghosts
    real, intent(inout) :: mass(n_real+n_ghosts), pos(n_real+n_ghosts), vel(n_real+n_ghosts), a(n_real+n_ghosts), &
    rho(n_real+n_ghosts), sml(n_real+n_ghosts), t,  u(n_real+n_ghosts), gamma
    real, intent(in) :: dt, dx
    real, intent(out) :: p(n_real+n_ghosts), cs(n_real+n_ghosts)
    integer, intent(inout) :: nfile
    real :: v_1(n_real + n_ghosts), a_0(n_real + n_ghosts)

    t = t+dt
    a_0(:) = a(:)


    pos(:) = pos(:) + dt*vel(:) + 0.5*(dt)**(2)*a(:)
    call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs, a)

    v_1(:) = vel(:) + dt*a(:)
    call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs, a)


    call derivs(n_real, n_ghosts, dx, pos, v_1, mass, rho, u, sml, gamma, p, cs, a, nfile)


    vel(:) = v_1(:) + 0.5*dt*(a(:) - a_0(:))
    call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs, a)


  end subroutine leap_frog_integrator

end module time_stepping
