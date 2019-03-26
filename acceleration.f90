module acceleration

  use spline_kernal, only: W
  !use viscosity, only: q, vsig
  implicit none


contains
  subroutine get_accel(n_real, n_ghosts, mass, vel, cs, pos, sml, rho, p, a)

    integer, intent(in) :: n_real, n_ghosts
    real, intent(in) :: pos(n_real+n_ghosts), mass(n_real+n_ghosts), sml(n_real+n_ghosts), &
    rho(n_real+n_ghosts), p(n_real+n_ghosts), vel(n_real+n_ghosts), cs(n_real+n_ghosts)
    real, intent(out) :: a(n_real+n_ghosts)
    real :: qa, qb, dx_a, dx_b, ma, mb, sml_i
    real :: w1a, w2a, w1b, w2b, pos_i, rho_i, dir
    real :: k, l
    integer :: i, j

    qa = 0.
    qb = 0.

    a(:) = 0.

    do i = 1,n_real
      sml_i = sml(i)
      pos_i = pos(i)
      rho_i = rho(i)

      do j = 1,n_real+n_ghosts

        dx_a = abs(pos_i - pos(j))/sml_i
        call w(dx_a, w1a, w2a)
        dx_b = abs(pos_i - pos(j))/sml(j)
        call w(dx_b, w1b, w2b)



        if (pos_i /= pos(j)) then
          dir = (pos_i-pos(j))/abs(pos_i - pos(j))
        else
          dir = 0.
        end if

        ! vsig_a = vsig(cs(i), vel(i) - vel(j), dir)
        ! qa = q(rho(i), vsig_a, vel(i) - vel(j), dir)
        !
        ! vsig_b = vsig(cs(j), vel(i) - vel(j), dir)
        ! qb = q(rho(j), vsig_b, vel(i) - vel(j), dir)


        k = (p(i) + qa)/rho_i**2 * dir * w2a/((sml_i)**2)
        l = (p(j) + qb)/rho(j)**2 * dir * w2b/((sml(j))**2)
        a(i) = a(i) - mass(j)*(k + l)
        !a(i) = a(i) + mass(j) * (((p(i) + qa)/rho_i**2) * dir * w2a/((sml_i)**2) +
        !((p(j) + qb)/rho(j)**2) * dir * w2b/((sml(j))**2))
      end do

    end do


  end subroutine get_accel

end module acceleration
