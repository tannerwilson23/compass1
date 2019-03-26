module setup

  implicit none

contains

  subroutine init_conditions_standing(n, pos, vel, mass, sml, rho, u, p, cs, spac)

    integer, intent(in) :: n
    real, intent(out) :: pos(n), vel(n), mass(n), sml(n), rho(n), u(n), p(n), cs(n), spac
    real :: xmin, xmax, rho0, h, amp
    real, parameter :: pi = 4.*atan(1.)
    integer :: i

    xmin = 0.0
    xmax = 1.0
    spac = (xmax - xmin)/(n)
    amp = 10.**(-4)

    rho0 = 1.0
    h = 1.2*spac

    rho(:) = rho0
    sml(:) = h
    mass(:) = rho0*(xmax-xmin)/(n)


    do i = 1,n
      pos(i) = (i-1)*spac
      vel(i) = amp*sin(2*pi*pos(i)/(xmax-xmin))
    enddo

  end subroutine init_conditions_standing


  subroutine init_conditions_sod(n, n_bound, pos, vel, mass, sml, rho, u, p, cs, spac)

    integer, intent(in) :: n, n_bound
    real, intent(out) :: pos(n), vel(n), mass(n), sml(n), rho(n), u(n), p(n), cs(n), spac
    real :: xmin, xmax, rho0, h, amp
    real, parameter :: pi = 4.*atan(1.)
    integer :: i

    xmin = -0.5
    xmax = 0.5


    do i = 1,n
      if (i == 1) then
        pos(i) = xmin
      else if (pos(i-1) < 0) then
        pos(i) = pos(i-1) + 0.001
      else
        pos(i) = pos(i-1) + 0.01
      end if
    enddo

    vel(1:n_bound) = 0.
    vel(n - n_bound:n) = 0.

  end subroutine init_conditions_sod

end module setup
