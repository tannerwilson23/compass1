module energy
  implicit none

contains
  subroutine kin_tot(n_real, m, v, T)
    integer, intent(in) :: n_real
    real, intent(in) :: m(n_real), v(n_real)
    real, intent(out) :: T
    integer :: i

    !calculate total kinetic energy of the system, usualy calculation.
    T = 0.
    do i = 1, n_real
      T = T + 0.5*m(i)*(v(i))**2
    enddo

  end subroutine kin_tot

  subroutine u_tot(n_real, m, u, E_int)
    integer, intent(in) :: n_real
    real, intent(in) :: m(n_real), u(n_real)
    real, intent(out) :: E_int
    integer :: i
    E_int = 0.
    do i = 1, n_real
      E_int = E_int + m(i) * u(i)
    enddo
  end subroutine u_tot

end module energy
