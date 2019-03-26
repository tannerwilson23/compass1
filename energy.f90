module energy
  implicit none

contains
  subroutine kin_tot(n_real, m, v, T)
    integer, intent(in) :: n_real
    real, intent(in) :: m(n_real), v(n_real)
    real, intent(out) :: T
    integer :: i

    T = 0
    do i = 1, n_real
      T = T = m(i)*(v(i))**2
    enddo
