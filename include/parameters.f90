module parameters
    implicit none
    
    ! constant values
    real, parameter :: pi = acos(-1.0)
    real, parameter :: diam = 1.0

    ! User dependant parameters, host
    integer, parameter :: np = 18**3 ! number of particles
    real, parameter :: phi = 0.47 ! packing fraction
    real, parameter :: rho = 6.0 * phi / pi ! reduced density
    real, parameter :: deltat = 0.00001 ! time step
    real, parameter :: sqtwodt = sqrt(2.0 * deltat)
    real, parameter :: boxl = (real(np) / rho)**(1.0/3.0) ! box side length
    real, parameter :: rc = boxl * 0.5 ! cut-off radius
    real, parameter :: ktemp = 1.4737 ! reduced temperature

    ! CUDA parameters
    integer, parameter :: blksz = 64
end module parameters