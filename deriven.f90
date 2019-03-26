module deriven

  use eos, only : equation_of_state
  use density, only: get_density
  use acceleration, only: get_accel
  use smoother, only: smoothed
  use ghosts, only: set_ghosts
  use output, only: write_output

  implicit none

  contains
    subroutine derivs(n_real, n_ghosts, dx, pos, vel, mass, rho, u, sml, gamma, p, cs, a, nfile)

      integer, intent(in) :: n_real, n_ghosts, nfile
      real, intent(inout) :: pos(n_real+n_ghosts), mass(n_real+n_ghosts), vel(n_real+n_ghosts), gamma, &
      u(n_real+n_ghosts)
      real, intent(inout) :: rho(n_real+n_ghosts), sml(n_real+n_ghosts)
      real, intent(out) :: cs(n_real+n_ghosts), a(n_real+n_ghosts), p(n_real+n_ghosts)
      real :: h_fac, dx
      integer :: j

      h_fac = 1.2

      do j = 1,4
        call get_density(n_real, n_ghosts, pos, mass, rho, sml)
        !call smoothed(n_real, h_fac, sml, mass, rho)
        call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs,a)
      enddo

      call equation_of_state(n_real, rho, gamma, p, cs)
      call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs,a)
      call get_accel(n_real, n_ghosts, mass, vel, cs, pos, sml, rho, p, a)
      call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs,a)


    end subroutine derivs
end module deriven
