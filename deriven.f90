module deriven

  use eos, only : equation_of_state
  use density, only: get_density
  use acceleration, only: get_accel
  use smoother, only: smoothed
  use ghosts, only: set_ghosts

  implicit none

  contains
    subroutine derivs(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, gamma, p, cs, u, du)

      integer, intent(in) :: n_real, n_ghosts, n_bound
      real, intent(inout) :: pos(n_real+n_ghosts), mass(n_real+n_ghosts), vel(n_real+n_ghosts), gamma, &
      u(n_real+n_ghosts), du(n_real + n_ghosts), rho(n_real+n_ghosts), sml(n_real+n_ghosts)
      real, intent(out) :: cs(n_real+n_ghosts), a(n_real+n_ghosts), p(n_real+n_ghosts)
      real :: h_fac
      integer :: j

      ! Set the smoothing factor, given in lecture
      h_fac = 1.2

      ! allowing the smoothing length to vary according to density calc, (7).
      !run three times, seems to converge.
      do j = 1,3
        call get_density(n_real, n_ghosts, pos, mass, rho, sml)
        call smoothed(n_real, h_fac, sml, mass, rho)
        call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)
      enddo


      !no tricks, finding p from eos and then using that to find the accel.
      call equation_of_state(n_real, rho, gamma, p, cs, u)
      call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)
      call get_accel(n_real, n_ghosts, mass, vel, cs, pos, sml, rho, p, a, du)
      call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)


    end subroutine derivs
end module deriven
