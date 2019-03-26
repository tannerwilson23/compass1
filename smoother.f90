module smoother
  implicit none

contains
  subroutine smoothed(n, h_fac, sml, mass, rho)

    integer, intent(in) :: n
    real, intent(in) :: h_fac, mass(n), rho(n)
    real, intent(inout) :: sml(n)

    sml(:) = h_fac*(mass(:)/rho(:))

  end subroutine smoothed

end module smoother
