program main

  use allocate, only: allocations, pos, vel, mass, sml, rho, u, p, cs, a
  use setup, only: init_conditions_standing, init_conditions_sod
  use ghosts, only: set_ghosts
  use output, only: write_output
  use deriven, only: derivs
  use density, only: get_density
  use time_stepping, only: leap_frog_integrator
  use time_step, only: step_size


  implicit none
  integer :: n, nfile, i, n_real, n_ghosts, choices, n_bound

  real :: gamma, dtout, time_cond, dt, tstop, t, dx

  print*, '1: Standing Wave', '2:Sod Shock 1', '3:Sod Shock 2?'
  read*, choices

  if (choices == 1) then
    print*, 'Hello World. How many cells would you like to use? '
    read*, n_real
  else if (choices == 2) then
    n_real = 551
  else
    n_real = 1000
  end if


  n_ghosts = 10
  n_ghosts = 2 * n_ghosts

  n = n_real + n_ghosts

  gamma = 1.
  dtout = 0.01
  t = 0.00001
  nfile = 0
  tstop = 5.
  time_cond = 0.2

  call allocations(n)

  if (choices == 1) then
    call init_conditions_standing(n_real, pos, vel, mass, sml, rho, u, p, cs, dx)
  else if (choices == 2) then
    n_bound = 6
    call init_conditions_sod(n_real, n_bound, pos, vel, mass, sml, rho, u, p, cs, dx)
  end if



  call write_output(nfile, n, pos, vel, mass, sml, rho, u, p, cs, a)

  call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs, a)
  call derivs(n_real, n_ghosts, dx, pos, vel, mass, rho, u, sml, gamma, p, cs, a, nfile)
  call write_output(nfile, n, pos, vel, mass, sml, rho, u, p, cs, a)

  call step_size(n_real, time_cond, sml, cs, dt)

   do i = 0,5000
    call leap_frog_integrator(n_real, n_ghosts, t, dx, dt, pos, vel, a, mass, rho, sml, gamma, p, cs, u, nfile)
    call set_ghosts(n_real, n_ghosts, dx, pos, vel, mass, sml, rho, u, p, cs, a)

    if (mod(t,dtout) < 1.e-2) then
      call write_output(nfile, n_real, pos, vel, mass, sml, rho, u, p, cs, a)
    end if

    call step_size(n_real, time_cond, sml, cs, dt)
  end do

end program main
