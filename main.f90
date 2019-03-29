program main

  ! Assignment 1 SPH solver
  ! Tanner Wilson
  ! Isothermal linear wave and Sod Shock Tube problem.

  ! import all useful routines from modules created
  use allocate, only: allocations, pos, vel, mass, sml, rho, u, p, cs, a, du
  use setup, only: init_conditions_standing, init_conditions_sod1, init_conditions_sod2
  use ghosts, only: set_ghosts
  use output, only: write_output
  use deriven, only: derivs
  use time_stepping, only: leap_frog_integrator, step_size


  implicit none

  ! Integers that are relatively self-explanatory
  integer :: n, nfile, n_real, n_ghosts, choices, n_bound, en_file_int, i

  ! Reals, gamma corresponds to adiabatic index
  ! dtout; how often we print output files
  ! tsto; when to stop program, depends on situation
  ! t; current time
  ! time_cond; condition on the time step related to courant condition
  real :: gamma, dtout, dtout_orig, time_cond, dt, tstop, t

  ! Which situation do we want to investigate
  print*, 'Hello World. Which Set-up? 1: Standing Wave   ', '2: Sod Shock 1   ', '3: Sod Shock 2?    '
  read*, choices

  ! Using the choice above, 1 has a choice of number of particles
  ! 2 and 3 do not.
  if (choices == 1) then
    print*, 'How many particles would you like to use? '
    read*, n_real
    ! Fixed number of ghosts, could change this to be another read but 10 each side
    ! is working well
    n_ghosts = 10
    n_ghosts = 2 * n_ghosts
    ! no boundary condition for the isothermal standing wave
    n_bound = 0
    gamma = 1.
    tstop = 6.
  else if (choices == 2) then
    ! fixed n_real particles
    n_real = 550
    n_ghosts = 10
    n_ghosts = 2 * n_ghosts
    ! number of boundary particles with v = 0
    n_bound = 6
    gamma = 1.0
    tstop = 0.3
  else
    n_real = 563
    gamma = 1.4
    tstop = 0.3
    n_ghosts = 10
    n_ghosts = 2 * n_ghosts
    ! number of boundary particles with v = 0
    n_bound = 6
  end if

  !total n
  n = n_real + n_ghosts

  !set various time variables
  dtout = 0.01
  dtout_orig = dtout
  t = 0.0

  time_cond = 0.2

  ! update nfile as the file number we are outputting
  nfile = 0

  !the energy unit, picked random large number that wont interfere with nfile
  en_file_int = 8000000

  ! allocate array sizes to be n
  call allocations(n)

  !self explanatory, set the initial conditions of the particles we are given
  if (choices == 1) then
    call init_conditions_standing(n_real, pos, vel, mass, sml)
  else if (choices == 2) then
    call init_conditions_sod1(n_real, pos, mass, sml)
  else
    call init_conditions_sod2(n_real, pos, mass, sml, u)
  end if

  ! open our energy file to write the energies to
  open(unit = en_file_int, file = 'outputs/energy', status = 'replace')
  write(en_file_int,*) 't   ', 'E_kin   ', 'E_int   ', 'E_tot   '

  ! First output of simply the input arguements, havent called derivs so rho etc = 0.
  call write_output(nfile, en_file_int, n, pos, vel, a, mass, rho, sml, p, cs, u, t)

  !set our ghost points so derivs works
  call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)


  !find our density, accel etc.
  call derivs(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, gamma, p, cs, u, du)

  ! write an output following finding rho, a, p etc. make sure there are no issues
  call write_output(nfile, en_file_int, n, pos, vel, a, mass, rho, sml, p, cs, u, t)

  ! find the minimum step size we can given the courant condition.
  call step_size(n_real, time_cond, sml, cs, dt)

  !Do a large number of times, i.e. much longer than the maximum time would correspond to.
   do i = 0,5000

     ! use our time stepping scheme
    call leap_frog_integrator(n_real, n_ghosts, n_bound, t, dt, pos, vel, a, mass, rho, sml, gamma, p, cs, u, du)

    !set the ghosts according to the new conditions; feel like im calling ghosts too many times.
    call set_ghosts(n_real, n_ghosts, n_bound, pos, vel, a, mass, rho, sml, p, cs, u, du)

    !essentially if the current time is a multiple of dtout then print an output

    if (t > dtout) then
      call write_output(nfile, en_file_int, n, pos, vel, a, mass, rho, sml, p, cs, u, t)
      dtout = dtout + dtout_orig
    end if

    !refind the step size given our new conditions
    call step_size(n_real, time_cond, sml, cs, dt)

    if (t > tstop) then
      stop
    end if
  end do


end program main
