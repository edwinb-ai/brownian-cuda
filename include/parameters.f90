module parameters
    use cudafor
    implicit none
    
    ! constant values
    real, parameter :: pi = acos(-1.0)
    real, parameter :: diam = 1.0

    ! User dependant parameters, host
    integer, parameter :: np = 18**3 ! number of particles
    real, parameter :: deltat = 0.00001 ! time step
    real, parameter :: sqtwodt = sqrt(2.0 * deltat)
    real, parameter :: ktemp = 1.4737 ! reduced temperature
    ! These variables can be read from a file
    real, managed :: boxl, rc
    real, managed :: phi, rho

    ! CUDA parameters
    integer, parameter :: blksz = 64
end module parameters