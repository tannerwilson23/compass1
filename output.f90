module output
  implicit none


  contains
    subroutine write_output(nfile ,n, pos, vel, mass, sml, rho, u, p, cs, a)

      integer, intent(in) :: n
      integer :: nfile
      real, intent(in) :: pos(n), vel(n), mass(n), sml(n), rho(n), u(n), p(n), cs(n), a(n)
      character(len = 40) :: filename
      integer :: i

      write(filename, "(a,i5.5)") 'outputs/output_', nfile

      open(nfile, file=filename, status = 'replace')

      do i = 1,n
        write(nfile,*) pos(i), vel(i), rho(i), p(i), cs(i), a(i), mass(i), sml(i)
      enddo
      close(nfile)
      nfile = nfile + 1

    end subroutine write_output

end module output
