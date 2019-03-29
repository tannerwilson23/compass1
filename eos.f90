module eos

  implicit none

contains
  subroutine equation_of_state(n, rho, gamma, p, cs, u)
    integer, intent(in) :: n
    real, intent(in) :: rho(n), gamma, u(n)
    real, intent(out) :: p(n), cs(n)

    !equation of state claculations, change for q17 to include u
    if (gamma == 1.) then
      p(:) = rho(:)**gamma
      cs(:) = sqrt(gamma * p(:)/rho(:))
    else
      p(:) = (gamma - 1.)*rho(:)*u(:)
      cs(:) = sqrt(gamma * p(:)/rho(:))
    end if

  end subroutine equation_of_state

end module eos
