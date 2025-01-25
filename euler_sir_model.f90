program lanier_seir_euler
    implicit none
    integer, parameter :: nruns = 10 ! number of runs
    integer, parameter :: npts = 10000 ! 10,000,000
    real(8), parameter :: days = 365.0d0 ! duration of simulation
    real(8) :: h, lam_sir, R0_sir, k_sir, lam_seir, R0_seir, k_seir, pop, S_sir(npts), I_sir(npts), R_sir(npts), E_seir(npts), I_seir(npts), R_seir(npts), time(npts) ! include arrays
    real(8) :: start1, finish1, start2, finish2, start3, finish3 ! timer
    real(8) :: total_program_time, total_method_time
    integer :: n, run
  
    call cpu_time(start3)
    ! initialize both timers
    total_program_time = 0.0d0
    total_method_time = 0.0d0
  
    do run = 1, nruns
      call cpu_time(start2) ! start program timer
  
      pop = 7900000.0d0 ! population
      R_sir(1) = 0.0d0 ! recovered individuals (SIR)
      I_sir(1) = 10.0d0 ! initial number of sick people (SIR)
      S_sir(1) = pop - I_sir(1) - R_sir(1)
      lam_sir = 7.0d0 ! recovery rate (duration of sickness in days) (SIR)
      R0_sir = 10.0d0 ! reproduction number (higher R0 means more contagious) NOTE: 1 person infects R0 more (SIR)
      k_sir = R0_sir / (S_sir(1) * I_sir(1) * 2.0d0) ! transmission rate (SIR)
  
      E_seir(1) = 0.0d0 ! initial number of exposed people (SEIR)
      I_seir(1) = 10.0d0 ! initial number of infected people (SEIR)
      R_seir(1) = 0.0d0 ! initial number of recovered individuals (SEIR)
      S_sir(1) = pop - E_seir(1) - I_seir(1) - R_seir(1)
      lam_seir = 5.0d0 ! rate of progression from exposed to infected (SEIR)
      R0_seir = 8.0d0 ! reproduction number (SEIR)
      k_seir = R0_seir / (S_sir(1) * I_seir(1) * 2.0d0) ! transmission rate (SEIR)
  
      h = days / npts ! timestep
  
      lam_sir = lam_sir * 2
      lam_seir = lam_seir * 2
  
      ! Euler (SEIR)
      call cpu_time(start1) ! start method timer
      do n = 2, npts
        S_sir(n) = S_sir(n - 1) + (-k_seir * S_sir(n - 1) * I_seir(n - 1)) * h
        E_seir(n) = E_seir(n - 1) + (k_seir * S_sir(n - 1) * I_seir(n - 1) - (1 / lam_seir) * E_seir(n - 1)) * h
        I_seir(n) = I_seir(n - 1) + (1 / lam_seir) * E_seir(n - 1) * h
        R_seir(n) = R_seir(n - 1) + (1 / lam_seir) * I_seir(n - 1) * h
      end do
  
      ! write results to files (SEIR)
      open(unit=400, file='s_seir.dat', status='replace') ! overwrites the .dat files
      open(unit=500, file='e_seir.dat', status='replace')
      open(unit=600, file='i_seir.dat', status='replace')
      open(unit=700, file='r_seir.dat', status='replace')
      ! headers
      write(400, *) "#t          s_seir"
      write(500, *) "#t          e_seir"
      write(600, *) "#t          i_seir"
      write(700, *) "#t          r_seir"
      ! write to files
      do n = 1, npts
        write(400, *) time(n), S_sir(n)
        write(500, *) time(n), E_seir(n)
        write(600, *) time(n), I_seir(n)
        write(700, *) time(n), R_seir(n)
      end do
      ! close files
      close(400)
      close(500)
      close(600)
      close(700)
  
      call cpu_time(finish1) ! end method timer
      total_method_time = (finish1 - start1)
  
      ! time vector
      do n = 1, npts
        time(n) = (n - 1) * h
      end do
  
      ! open files (SIR)
      open(unit=100, file='s_sir.dat', status='replace') ! overwrites the .dat files
      open(unit=200, file='i_sir.dat', status='replace')
      open(unit=300, file='r_sir.dat', status='replace')
      ! headers
      write(100, *) "#t          s_sir"
      write(200, *) "#t          i_sir"
      write(300, *) "#t          r_sir"
      ! write to files (SIR)
      do n = 1, npts
        write(100, *) time(n), S_sir(n)
        write(200, *) time(n), I_sir(n)
        write(300, *) time(n), R_sir(n)
      end do
      ! close files (SIR)
      close(100)
      close(200)
      close(300)
  
      call cpu_time(finish2) ! end program timer
      total_program_time = (finish2 - start2)
      ! write(*,*) 'run: ', run ,'program runtime', (total_program_time)
    end do
  
    call cpu_time(finish3)
    write(*,*) 'Duration of simulation:', days, 'days.'
    write(*,*) 'Number of points:', npts
    write(*,*) 'Total population:', pop, 'people.'
    write(*,*) 'Initial sick (SIR):', I_sir(1), 'people.'
    write(*,*) 'Initial infected (SEIR):', I_seir(1), 'people.'
  
    write(*,*) 'Average method runtime: ', (total_method_time) / real(nruns), 'seconds'
    write(*,*) 'Average program runtime: ', (total_program_time) / real(nruns), 'seconds'
    write(*,*) 'Total program runtime: ', (finish3 - start3), 'seconds'
  end program lanier_seir_euler
  