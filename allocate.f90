module allocate

  implicit none

  real, allocatable :: pos(:), vel(:), mass(:), sml(:), rho(:), p(:), cs(:), a(:), u(:), du(:)


  public :: pos, vel, mass, sml, rho, u, p, cs, du

contains
  subroutine allocations(n)
    integer, intent(in) :: n

    !allocate arrays to be of size n
    allocate(pos(n), vel(n), mass(n), sml(n), rho(n), p(n), cs(n), a(n), u(n), du(n))

  end subroutine allocations

end module allocate
