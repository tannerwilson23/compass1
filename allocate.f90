module allocate

  implicit none

  real, allocatable :: pos(:), vel(:), mass(:), sml(:), rho(:), u(:), p(:), cs(:), a(:)


  public :: pos, vel, mass, sml, rho, u, p, cs

contains
  subroutine allocations(n)
    integer, intent(in) :: n

    allocate(pos(n), vel(n), mass(n), sml(n), rho(n), u(n), p(n), cs(n), a(n))

  end subroutine allocations

end module allocate
