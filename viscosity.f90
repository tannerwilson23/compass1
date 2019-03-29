module viscosity
  implicit none
  public :: vsig, q

  !artifical viscosity calculations according to (8) and (9) in q12.

contains
    real function vsig(csa, vab, rab)
      implicit none
      real, intent(in) :: csa, vab, rab
      real :: a, b
      a = 0
      b = 0

      vsig = a*csa - b*(vab*rab)

    end function vsig

    real function q(rhoa, vsiga, vab, rab)
      implicit none
      real, intent(in) :: rhoa, vsiga, vab, rab

      if (vab * rab < 0) then
        q = -0.5*rhoa*vsiga*vab*rab
      else
        q = 0.
      end if


    end function q

end module viscosity
