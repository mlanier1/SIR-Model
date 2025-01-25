program lanier_sir_rk2
  implicit none
  integer, parameter :: nruns = 10 ! number of runs
  integer, parameter :: npts = 10000 !10000000
  real(8), parameter :: days = 365.0d0 ! duration of simulation
  real(8) :: h, lam, R0, k, pop, S(npts), I(npts), R(npts), time(npts) ! include arrays
  real(8) :: dS1, dI1, dR1, S_mid, I_mid, R_mid, dS2, dI2, dR2 ! derivatives/midpoints
  real(8) :: start1, finish1, start2, finish2, start3, finish3 ! timer
  real(8) :: total_program_time, total_method_time
  integer :: n, run

  call cpu_time(start3)
  ! initialize both timers
  total_program_time = 0.0d0
  total_method_time = 0.0d0

  !do run = 1, nruns
    call cpu_time(start1) ! start program timer
    ! write(*,*) 'Duration of simulation:',days,'days.'
    ! write(*,*) 'Number of points:',npts

    pop = 7900000.0d0 ! population
    R(1) = 0.0d0 ! recovered individuals
    I(1) = 10.0d0 ! start with one sick person
    S(1) = pop - I(1) - R(1)
    lam = 7.0d0 ! recovery rate (duration of sickness in days)
    R0 = 10.0d0 ! reproduction number (higher R0 means more contagious)
    k = R0 / (S(1) * I(1) * 2.0d0) ! transmission rate
    h = days / npts ! timestep
lam = lam * 2
    ! write(*,*) 'Total population:',pop,'people.'
    ! write(*,*) 'Initial sick:',I(1),'people.'
    do run = 1, nruns

      ! RK2
      call cpu_time(start2) ! start method timer
      do n = 2, npts
        ! first set of derivatives
        dS1 = -k * S(n - 1) * I(n - 1)
        dI1 = k * S(n - 1) * I(n - 1) - (1.0d0 / lam) * I(n - 1)
        dR1 = (1.0d0 / lam) * I(n - 1)
        ! derivatives to estimate values at the midpoint
        S_mid = S(n - 1) + dS1 * h / 2.0d0
        I_mid = I(n - 1) + dI1 * h / 2.0d0
        R_mid = R(n - 1) + dR1 * h / 2.0d0
        ! second set of derivatives using the estimated values
        dS2 = -k * S_mid * I_mid
        dI2 = k * S_mid * I_mid - (1.0d0 / lam) * I_mid
        dR2 = (1.0d0 / lam) * I_mid
        ! second set of derivatives to update the variables
        S(n) = S(n - 1) + dS2 * h
        I(n) = I(n - 1) + dI2 * h
        R(n) = R(n - 1) + dR2 * h
      end do
      call cpu_time(finish2)
      ! write(*,'(A,ES10.3,A)') ' Method runtime: ', (finish2 - start2), 'seconds.'
      total_method_time = (finish2-start2)

      ! time vector
      do n = 1, npts
        time(n) = (n - 1) * h
      end do

      ! open files
    open(unit=100, file='s_sir_rk2.dat', status='replace') ! overwrites the .dat files
    open(unit=200, file='i_sir_rk2.dat', status='replace')
    open(unit=300, file='r_sir_rk2.dat', status='replace')
    ! headers
    write(100, *) "#t          s"
    write(200, *) "#t          i"
    write(300, *) "#t          r"
    ! write to files
    do n = 1, npts
      write(100, *) time(n), S(n)
      write(200, *) time(n), I(n)
      write(300, *) time(n), R(n)
    end do
    ! close files
    close(100)
    close(200)
    close(300)

  call cpu_time(finish1) ! end program timer
  total_program_time = (finish1-start1)
end do
  call cpu_time(finish3)
  write(*,*) 'Duration of simulation:',days,'days.'
  write(*,*) 'Number of points:',npts
  write(*,*) 'Total population:',pop,'people.'
  write(*,*) 'Initial sick:',I(1),'people.'
  
  write(*,*) 'Average method runtime: ', (total_method_time)/real(nruns), 'seconds'
  write(*,*) 'Average program runtime: ', (total_program_time)/real(nruns), 'seconds'
  write(*,*) 'Total program runtime: ', (finish3-start3)
end program lanier_sir_rk2
