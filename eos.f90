module eos

  implicit none

contains
  subroutine equation_of_state(n, rho, gamma, p, cs)
    integer, intent(in) :: n
    real, intent(in) :: rho(n), gamma
    real, intent(out) :: p(n), cs(n)

    p(:) = rho(:)**gamma
    cs(:) = sqrt(gamma * p(:)/rho(:))

  end subroutine equation_of_state

end module eos
