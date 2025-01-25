program lanier_sir_rk4
  implicit none
  integer, parameter :: nruns = 10 ! number of runs
  integer, parameter :: npts = 10000
  real(8), parameter :: days = 365.0d0 ! duration of simulation
!  real(8), parameter :: pop = 100000.0d0 ! population
 ! real(8), parameter :: lam = 8.0d0 ! recovery rate (duration of sickness in days)
  !real(8), parameter :: R0 = 10.0d0 ! reproduction number (higher R0 means more contagious)
  real(8) :: h, lam, k, R0,pop, S(npts), I(npts), R(npts), time(npts) ! include arrays
  real(8) :: S1, I1, R1, S2, I2, R2, S3, I3, R3, S4, I4, R4 ! each stage
  real(8) :: dS1, dI1, dR1, dS2, dI2, dR2, dS3, dI3, dR3, dS4, dI4, dR4 ! derivatives/slopes
  real(8) :: start1, finish1, start2, finish2, start3, finish3 ! timer
  real(8) :: total_program_time, total_method_time
  integer :: n, run

call cpu_time(start3)
! initialize both timers
total_program_time = 0.0d0
total_method_time = 0.0d0

call cpu_time(start1) ! start program timer
!  write(*,*) 'Duration of simulation:', days, 'days.'
 ! write(*,*) 'Number of points:', npts

  pop = 7900000.0d0 ! population
  R(1) = 0.0d0 ! recovered individuals
  I(1) = 10.0d0 ! start with one sick person
  S(1) = pop - I(1) - R(1)
  lam = 7.0d0 ! recovery rate (duration of sickness in days)
  R0 = 10.0d0 ! reproduction number (higher R0 means more contagious)
  k = R0 / (S(1) * I(1) * 2.0d0) ! transmission rate
  h = days / npts ! timestep
lam = lam * 2
!  write(*,*) 'Total population:', pop, 'people.'
 ! write(*,*) 'Initial sick:', I(1), 'people.'
do run = 1, nruns
  ! RK4
  call cpu_time(start2)
  do n = 2, npts

    dS1 = -k * S(n - 1) * I(n - 1)
    dI1 = k * S(n - 1) * I(n - 1) - (1.0d0 / lam) * I(n - 1)
    dR1 = (1.0d0 / lam) * I(n - 1)
    
    S1 = S(n - 1) + dS1 * h / 2.0d0
    I1 = I(n - 1) + dI1 * h / 2.0d0
    R1 = R(n - 1) + dR1 * h / 2.0d0

    dS2 = -k * S1 * I1
    dI2 = k * S1 * I1 - (1.0d0 / lam) * I1
    dR2 = (1.0d0 / lam) * I1

    S2 = S(n - 1) + dS2 * h / 2.0d0
    I2 = I(n - 1) + dI2 * h / 2.0d0
    R2 = R(n - 1) + dR2 * h / 2.0d0

    dS3 = -k * S2 * I2
    dI3 = k * S2 * I2 - (1.0d0 / lam) * I2
    dR3 = (1.0d0 / lam) * I2

    S3 = S(n - 1) + dS3 * h
    I3 = I(n - 1) + dI3 * h
    R3 = R(n - 1) + dR3 * h

    dS4 = -k * S3 * I3
    dI4 = k * S3 * I3 - (1.0d0 / lam) * I3
    dR4 = (1.0d0 / lam) * I3

    ! update variables
    S(n) = S(n - 1) + (dS1 + 2.0d0 * dS2 + 2.0d0 * dS3 + dS4) * h / 6.0d0
    I(n) = I(n - 1) + (dI1 + 2.0d0 * dI2 + 2.0d0 * dI3 + dI4) * h / 6.0d0
    R(n) = R(n - 1) + (dR1 + 2.0d0 * dR2 + 2.0d0 * dR3 + dR4) * h / 6.0d0
  end do

  call cpu_time(finish2)
  total_method_time = (finish2-start2)
  !write(*,'(A,ES10.3,A)') ' Method runtime: ', (finish2 - start2), 'seconds.'

  ! time vector
  do n = 1, npts
    time(n) = (n - 1) * h
  end do

  ! open files
  open(unit=100, file='s_sir_rk4.dat', status='replace') ! overwrites the .dat files
  open(unit=200, file='i_sir_rk4.dat', status='replace')
  open(unit=300, file='r_sir_rk4.dat', status='replace')

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
end program lanier_sir_rk4
