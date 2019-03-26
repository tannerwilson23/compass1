module ghosts
  implicit none

contains

  subroutine set_ghosts(n_real, g_n, dx, pos, vel, mass, sml, rho, u, p, cs, a)
    integer, intent(in) :: n_real, g_n
    real, intent(in) :: dx
    real, intent(inout) :: pos(n_real+g_n), vel(n_real+g_n), mass(n_real+g_n), &
    sml(n_real+g_n), rho(n_real+g_n), u(n_real+g_n), p(n_real+g_n), cs(n_real+g_n), &
    a(n_real+g_n)

    if (pos(n_real + 1) == 0) then
      pos(n_real+g_n/2+1:n_real+g_n) = pos(1:1+g_n/2)-(g_n/2)*dx
      pos(n_real+1:n_real+g_n/2) = pos(n_real-g_n/2+1:n_real)+((g_n/2)*dx)
    end if

    vel(n_real+1:n_real+g_n/2) = vel(1:g_n/2)
    mass(n_real+1:n_real+g_n/2) = mass(1:g_n/2)
    sml(n_real+1:n_real+g_n/2) = sml(1:g_n/2)
    rho(n_real+1:n_real+g_n/2) = rho(1:g_n/2)
    u(n_real+1:n_real+g_n/2) = u(1:g_n/2)
    p(n_real+1:n_real+g_n/2) = p(1:g_n/2)
    cs(n_real+1:n_real+g_n/2) = cs(1:g_n/2)
    a(n_real+1:n_real+g_n/2) = a(1:g_n/2)





    vel(n_real+g_n/2+1:n_real+g_n) = vel(n_real-g_n/2+1:n_real)
    mass(n_real+g_n/2+1:n_real+g_n) = mass(n_real-g_n/2+1:n_real)
    sml(n_real+g_n/2+1:n_real+g_n) = sml(n_real-g_n/2+1:n_real)
    rho(n_real+g_n/2+1:n_real+g_n) = rho(n_real-g_n/2+1:n_real)
    u(n_real+g_n/2+1:n_real+g_n) = u(n_real-g_n/2+1:n_real)
    p(n_real+g_n/2+1:n_real+g_n) = p(n_real-g_n/2+1:n_real)
    a(n_real+g_n/2+1:n_real+g_n) = a(n_real-g_n/2+1:n_real)

  end subroutine set_ghosts



end module ghosts
