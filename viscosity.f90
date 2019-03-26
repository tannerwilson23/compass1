module viscosity
  implicit none

  contains

    function vsig(csa, vab, rab)
      real, intent(in) :: csa, vab, rab
      real :: a, b
      real, intent(out) :: vsig

      vsig = a*csa - b(vab*rab)

      return vsig

    function q(rhoa, vsiga, vab, rab)
      real, intent(in) :: rhoa, vsiga, vab, rab
      real, intent(out) :: q

      if (vab * rab < 0) then
        q = -0.5*rhoa*vsiga*vab*rab
      else
        q = 0.
      
      return q
