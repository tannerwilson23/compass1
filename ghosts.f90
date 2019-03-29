module ghosts
  implicit none

contains

  subroutine set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)
    integer, intent(in) :: n_real, n_ghosts, n_bound
    real, intent(inout) :: pos(n_real+n_ghosts), vel(n_real+n_ghosts), mass(n_real+n_ghosts), &
    sml(n_real+n_ghosts), rho(n_real+n_ghosts), p(n_real+n_ghosts), cs(n_real+n_ghosts), &
    a(n_real+n_ghosts), u(n_real+n_ghosts), du(n_real+n_ghosts)

    !the first n_ghost/2 correspond to those on the left of the domain and the second on the right
    !n_ghost set to n_ghost*2 in main to avoid halving integers.

    !if the positions of the ghost particles haven't been set, set them.
    !updating the positions every time was causing clumping of particles, density issues
    ! shocks. let time stepping update ghost points.
    if (pos(n_real + 1) == 0) then
      if (n_real == 550) then
        !if the sod shock problem, continue the trend of the spacing depending on the side you are on
        pos(n_real+1:n_real+n_ghosts/2) = pos(1:n_ghosts/2) - n_ghosts/2 * 0.001
        pos(n_real+n_ghosts/2+1:n_real+n_ghosts) = pos(n_real-n_ghosts/2+1:n_real) + n_ghosts/2 * 0.01
      else if (n_real == 563) then
        pos(n_real+1:n_real+n_ghosts/2) = pos(1:n_ghosts/2) - n_ghosts/2 * 0.001
        pos(n_real+n_ghosts/2+1:n_real+n_ghosts) = pos(n_real-n_ghosts/2+1:n_real) + n_ghosts/2 * 0.008
      else
        pos(n_real+n_ghosts/2+1:n_real+n_ghosts) = pos(n_real-n_ghosts/2+1:n_real) - 1.
        pos(n_real+1:n_real+n_ghosts/2) = pos(1:1+n_ghosts/2) + 1.
      end if
    end if

    !probably a more elegant way of doing these but setting the ghost particle conditions
    !according to the first n_ghosts/2 or n_ghosts particles

    vel(n_real+1:n_real+n_ghosts/2) = vel(1:n_ghosts/2)
    mass(n_real+1:n_real+n_ghosts/2) = mass(1:n_ghosts/2)
    sml(n_real+1:n_real+n_ghosts/2) = sml(1:n_ghosts/2)
    rho(n_real+1:n_real+n_ghosts/2) = rho(1:n_ghosts/2)
    u(n_real+1:n_real+n_ghosts/2) = u(1:n_ghosts/2)
    p(n_real+1:n_real+n_ghosts/2) = p(1:n_ghosts/2)
    cs(n_real+1:n_real+n_ghosts/2) = cs(1:n_ghosts/2)
    a(n_real+1:n_real+n_ghosts/2) = a(1:n_ghosts/2)
    du(n_real+1:n_real+n_ghosts/2) = du(1:n_ghosts/2)

    vel(n_real+n_ghosts/2+1:n_real+n_ghosts) = vel(n_real-n_ghosts/2+1:n_real)
    mass(n_real+n_ghosts/2+1:n_real+n_ghosts) = mass(n_real-n_ghosts/2+1:n_real)
    sml(n_real+n_ghosts/2+1:n_real+n_ghosts) = sml(n_real-n_ghosts/2+1:n_real)
    rho(n_real+n_ghosts/2+1:n_real+n_ghosts) = rho(n_real-n_ghosts/2+1:n_real)
    cs(n_real+n_ghosts/2+1:n_real+n_ghosts) = cs(n_real-n_ghosts/2+1:n_real)
    u(n_real+n_ghosts/2+1:n_real+n_ghosts) = u(n_real-n_ghosts/2+1:n_real)
    p(n_real+n_ghosts/2+1:n_real+n_ghosts) = p(n_real-n_ghosts/2+1:n_real)
    a(n_real+n_ghosts/2+1:n_real+n_ghosts) = a(n_real-n_ghosts/2+1:n_real)
    du(n_real+n_ghosts/2+1:n_real+n_ghosts) = du(n_real-n_ghosts/2+1:n_real)

    !set the boundary conditions in the sod shock tube problem, could do this in another
    !subroutines though I feel like I would call it as much as set_ghosts
    if (n_bound /= 0) then
      vel(1:n_bound) = 0.
      vel(n_real-n_bound:n_real) = 0.
    end if

  end subroutine set_ghosts



end module ghosts
