program main
    use parameters
    use cudafor
    use curand
    use randomf
    use forcec, only: force
    use positions, only: position
    use utils
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    ! CUDA variables
    type(dim3) :: grid, tBlock
    integer :: numSMs
    integer :: istat
    type(curandGenerator) :: gen
    
    ! General arrays, either CUDA or host
    real, allocatable, device :: randv(:)
    real, allocatable, device :: fx_d(:), fy_d(:), fz_d(:)
    real, allocatable, managed :: x(:), y(:), z(:)
    real, allocatable, managed :: enerpot(:), zfac(:)

    ! Local variables
    real :: d, virial, epotn
    integer :: i, istep, nprom, j, u
    real :: dv, fnorm, pressure, sqtwodt, compressibility
    integer :: limT, limG, pbc, avefreq, threenp

    ! Variable initialization
    threenp = 3 * np
    d = (1.0 / rho)**(1.0/3.0) ! Interparticle distance
    avefreq = 10000 ! Average frequency
    sqtwodt = sqrt(2.0 * deltat)
    ! Always start with PBC
    pbc = 1
    ! Thermalization and equilibrium steps
    limT = 1500000
    limG = 3000000

    ! Initialize arrays
    allocate(x(np), y(np), z(np), source=0.0)
    allocate(enerpot(np), zfac(np), source=0.0)
    allocate(randv(threenp), source=0.0)
    allocate(fx_d(np), fy_d(np), fz_d(np), source=0.0)

    ! Show important information to screen
    write(unit=output_unit, fmt='(a, 8f)'), 'The length of the box is: ', boxl
    write(unit=output_unit, fmt='(a, 8f)'), 'The mean interparticle distance is: ', d
    write(unit=output_unit, fmt='(a, 8f)'), 'Cut radius: ', rc
    write(unit=output_unit, fmt='(a, 8f)'), 'Reduced temperature: ', ktemp

    ! Initialize CUDA kernel variables
    tBlock = dim3(blksz, 1, 1)
    istat = cudaDeviceGetAttribute(numSMs, cudaDevAttrMultiProcessorCount, 0)
    grid = 32 * numSMs
    write(unit=output_unit, fmt='(2i8)') grid%x, numSMs
    call force <<< grid, tBlock >>> (x, y, z, fx_d, fy_d, fz_d, enerpot, zfac)
    istat = cudaDeviceSynchronize()
    if ( istat .ne. 0 ) then
        write(unit=u, fmt='(a)') 'Error with GPU!'
    end if

    ! Create a new configuration
    call iniconfig(x, y, z, d)
    ! Energy of the initial configuration
    write(unit=output_unit, fmt='(a, 8f)'), 'E/N=', sum(enerpot)/real(np)

    ! Initialize the PRNG
    call initialize_rng(gen)

    open(newunit=u, file='energy_BD_therm.dat', status='unknown')

    write(unit=output_unit, fmt='(a)') 'Step, Upot, Press, Z'
    do istep = 1, limT
        ! Compute the random numbers
        call fill_vec(randv, gen, threenp, sqtwodt)

        call position <<< grid, tBlock >>> (x, y, z, fx_d, fy_d, fz_d, randv, pbc)
        call force <<< grid, tBlock >>> (x, y, z, fx_d, fy_d, fz_d, enerpot, zfac)

        if (mod(istep, avefreq) == 0) then
            ! Compute the energy per particle
            epotn = sum(enerpot) / real(np)
            
            ! Compute the virial, 3.0 is because of the dimensionality
            virial = sum(zfac) / 3.0
            
            ! Compute the pressure from the system
            pressure = (rho * ktemp) + (virial / boxl**3)

            ! Compute the compressibility factor
            compressibility = pressure / (rho * ktemp)

            ! Print to screen
            write(unit=output_unit, fmt='(i,4f16.8)') istep, epotn, pressure, &
                compressibility
            ! Print to file 'u'
            write(unit=u, fmt='(i,3f16.8)') istep, epotn, pressure, compressibility
        end if
    end do
    close (u)

    write(unit=u, fmt='(a)'), 'The system has thermalized'

    open(newunit=u, file='finalconBD.dat', status='unknown')
    do i = 1, np
        write(unit=u, '(3f16.8)') x(i), y(i), z(i)
    end do
    close(u)

    ! Open file for production run
    open(newunit=u, file='energy_BD_prod.dat', status='unknown')
    
    write(unit=output_unit, fmt='(a)') 'Step, Upot, Press, Z'
    do istep = 1, limG
        ! Compute the random numbers
        call fill_vec(randv, gen, threenp, sqtwodt)

        call position <<< grid, tBlock >>> (x, y, z, fx_d, fy_d, fz_d, randv, pbc)
        call force <<< grid, tBlock >>> (x, y, z, fx_d, fy_d, fz_d, enerpot, zfac)

        if (mod(istep, 1000) == 0) then
            ! Compute the energy per particle
            epotn = sum(enerpot) / real(np)
            
            ! Compute the virial, 3.0 is because of the dimensionality
            virial = sum(zfac) / 3.0
            
            ! Compute the pressure from the system
            pressure = (rho * ktemp) + (virial / boxl**3)

            ! Compute the compressibility factor
            compressibility = pressure / (rho * ktemp)

            ! Print to screen
            write(unit=output_unit, fmt='(i,4f16.8)') istep, epotn, pressure, &
                compressibility
            ! Print to file 'u'
            write(unit=u, fmt='(i,3f16.8)') istep, epotn, pressure, compressibility
        end if
    end do
    
    ! Close file for production
    close(u)
end program main
