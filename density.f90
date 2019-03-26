module density

use spline_kernal, only: W

implicit none

contains
  subroutine get_density(n_real, n_ghosts, pos, mass, rho, sml)

    integer, intent(in) :: n_real, n_ghosts
    real, intent(in) :: pos(n_real+n_ghosts), mass(n_real+n_ghosts), sml(n_real+n_ghosts)
    real, intent(out) :: rho(n_real+n_ghosts)
    real :: q, w_val, w_grad
    integer :: i,j


    rho(:) = 0

    do i = 1,n_real
      do j = 1,n_real+n_ghosts
        q = abs(pos(j) - pos(i))/sml(i)
        
        call w(q, w_val, w_grad)
        rho(i) = rho(i) + mass(j)*w_val/sml(i)
      enddo
    enddo

  end subroutine get_density


end module density
