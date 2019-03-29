module acceleration

  use spline_kernal, only: W
  use viscosity, only: q, vsig
  implicit none


contains
  subroutine get_accel(n_real, n_ghosts, mass, vel, cs, pos, sml, rho, p, a, du)

    integer, intent(in) :: n_real, n_ghosts
    real, intent(in) :: pos(n_real+n_ghosts), mass(n_real+n_ghosts), sml(n_real+n_ghosts), &
    rho(n_real+n_ghosts), p(n_real+n_ghosts), vel(n_real+n_ghosts), cs(n_real+n_ghosts)
    real, intent(out) :: a(n_real+n_ghosts), du(n_real+n_ghosts)
    real :: qa, qb, dx_a, dx_b, sml_i
    real :: w1a, w2a, w1b, w2b, pos_i, rho_i, dir, vsig_a, vsig_b
    real :: k, l, velij
    integer :: i, j

    !set qa and qb. to be zero, so if the viscosity is turned on they will be updated
    ! if not they could take on random values which will mess with this
    qa = 0.
    qb = 0.

    !set accl and u "accel" to be zero before summing
    a(:) = 0.
    du(:) = 0.

    !update only our real particle values, use set ghosts to set others
    do i = 1, n_real
      sml_i = sml(i)
      pos_i = pos(i)
      rho_i = rho(i)

      !sum over ghost points and real points for their contriubtion
      do j = 1,n_real+n_ghosts

        !find the gradW values for the two q values, (q in paper not viscosity)
        dx_a = abs(pos_i - pos(j))/sml_i
        call w(dx_a, w1a, w2a)
        dx_b = abs(pos_i - pos(j))/sml(j)
        call w(dx_b, w1b, w2b)

        !find the unit vector of the difference, so we don't get 0/0 set to 0 if equal
        if (pos_i /= pos(j)) then
          dir = (pos_i-pos(j))/abs(pos_i - pos(j))
        else
          dir = 0.
        end if

        !velocity difference, used a few times so just set variable.
        velij = vel(i) - vel(j)

        !aritificial viscosity from Q12, functions from viscosity module
        vsig_a = vsig(cs(i), velij, dir)
        qa = q(rho(i), vsig_a, velij, dir)

        vsig_b = vsig(cs(j), velij, dir)
        qb = q(rho(j), vsig_b, velij, dir)

        !k and l just the two different components in the calcualtion, divide and conquer, don't do in one step
        k = (p(i) + qa)/rho_i**2 * dir * w2a/((sml_i)**2)
        l = (p(j) + qb)/rho(j)**2 * dir * w2b/((sml(j))**2)
        a(i) = a(i) - mass(j)*(k + l)


        !update du according to (10)
        du(i) = du(i) + mass(j) * (p(i) + qa)/(rho_i**2) * (velij) * dir *  w2a /sml_i**2

      end do

    end do


  end subroutine get_accel

end module acceleration
