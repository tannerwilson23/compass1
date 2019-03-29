module setup

  implicit none

contains

  !conditions according to q2.
  subroutine init_conditions_standing(n, pos, vel, mass, sml)

    integer, intent(in) :: n
    real, intent(out) :: pos(n), vel(n), mass(n), sml(n)
    real :: xmin, xmax, rho0, h, amp, spac
    real, parameter :: pi = 4.*atan(1.)
    integer :: i

    xmin = 0.0
    xmax = 1.0
    spac = (xmax - xmin)/(n)

    amp = 10.**(-4) !* cs = 1. in the isothermal case

    h = 1.2*spac
    sml(:) = h

    !mass is based upon the density
    rho0 = 1.0
    mass(:) = rho0*(xmax-xmin)/(n)

    !equally spaced particles, velocity is a sin wave
    do i = 1,n
      pos(i) = (i-1)*spac
      vel(i) = amp*sin(2*pi*pos(i)/(xmax-xmin))
    enddo

  end subroutine init_conditions_standing

  !sod shock tube set up
  subroutine init_conditions_sod1(n, pos, mass, sml)

    integer, intent(in) :: n
    real, intent(out) :: pos(n), mass(n), sml(n)
    real :: xmin, xmax
    real, parameter :: pi = 4.*atan(1.)
    integer :: i

    xmin = -0.5
    xmax = 0.5

    mass(:) = 0.5*1.0/500

    !pretty self explanatory; if pos(i-1) < 0 then continue to add particles with
    !small spacing until you reach 0 and then make more spread out.
    do i = 1,n
      if (i == 1) then
        pos(i) = xmin
        sml(i) = 1.2 * 0.001
      else if (pos(i-1) < 0) then
        pos(i) = pos(i-1) + 0.001
        sml(i) = 1.2 * 0.001
      else
        pos(i) = pos(i-1) + 0.01
        sml(i) = 1.2 * 0.01
      end if
    enddo


  end subroutine init_conditions_sod1

  subroutine init_conditions_sod2(n, pos, mass, sml, u)

    integer, intent(in) :: n
    real, intent(out) :: pos(n), mass(n), sml(n), u(n)
    real :: xmin, xmax
    integer :: i

    xmin = -0.5
    xmax = 0.5

    mass(:) = 0.5*1.0/500

    !pretty self explanatory; if pos(i-1) < 0 then continue to add particles with
    !small spacing until you reach 0 and then make more spread out.
    do i = 1,n
      if (i == 1) then
        pos(i) = xmin
        sml(i) = 1.2 * 0.001
        u(i) = 2.5
      else if (pos(i-1) < 0) then
        pos(i) = pos(i-1) + 0.001
        sml(i) = 1.2 * 0.001
        u(i) = 2.5
      else
        pos(i) = pos(i-1) + 0.008
        sml(i) = 1.2 * 0.008
        u(i) = 2.
      end if
    enddo


  end subroutine init_conditions_sod2

end module setup
