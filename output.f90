module output

  use energy, only: kin_tot, u_tot
  implicit none


  contains
    subroutine write_output(nfile, en_file_int, n, pos, vel, a, mass, rho, sml, p, cs, u, t)

      integer, intent(in) :: n, en_file_int
      integer, intent(inout) :: nfile
      real, intent(in) :: pos(n), vel(n), mass(n), sml(n), rho(n), u(n), p(n), cs(n), a(n), t
      real :: E_kin, E_int
      character(len = 40) :: filename
      integer :: i

      !write the filename to a string
      write(filename, "(a,i5.5)") 'outputs/output_', nfile

      !open the file
      open(nfile, file=filename, status = 'replace')

      !appropriate column names
      write(nfile,*) 'x ', 'v ', '/rho  ', 'P ', 'cs  ', 'a ', 'mass  ', 'h '
      write(nfile,*) t

      !for all of our particles, write the value of the variable for that given particle
      do i = 1,n
        write(nfile,*) pos(i), vel(i), rho(i), p(i), cs(i), a(i), mass(i), sml(i), u(i)
      enddo
      close(nfile)

      !update the file name
      nfile = nfile + 1

      !write the kinetic energy to the energy file
      call u_tot(n, mass, u, E_int)
      call kin_tot(n, mass, vel, E_kin)
      write(en_file_int,*) t, E_kin, E_int, E_kin + E_int

    end subroutine write_output



end module output
