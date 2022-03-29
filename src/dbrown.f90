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
    real :: dv, fnorm, pressure, compressibility
    integer :: limT, limG, pbc, avefreq, threenp
    logical :: exists, init_config

    !! Read the inputs from the file
    call parse_input("input.in", limG, limT, init_config)
    ! Variable initialization
    rho = 6.0 * phi / pi
    boxl = (real(np) / rho)**(1.0/3.0)
    rc = boxl / 2.0
    threenp = 3 * np
    d = (1.0 / rho)**(1.0/3.0) ! Interparticle distance
    avefreq = 10000 ! Average frequency
    ! Always start with PBC
    pbc = 1

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
    write(unit=output_unit, fmt='(a, i)'), 'Number of particles: ', np

    ! Initialize CUDA kernel variables
    tBlock = dim3(blksz, 1, 1)
    istat = cudaDeviceGetAttribute(numSMs, cudaDevAttrMultiProcessorCount, 0)
    grid = 32 * numSMs
    write(unit=output_unit, fmt='(2i8)') grid%x, numSMs

    ! Check if there is a previous configuration
    inquire(file='finalconBD.dat', exist=exists)
    ! If so, load it
    if (exists .and. init_config) then
        open(newunit=u, file='finalconBD.dat', status="old", action="read")
        write(unit=output_unit, fmt='(a)') 'Reading positions from file...'
        do i = 1, np
            read(u, *) x(i), y(i), z(i)
        end do
        close(u)
    ! If not, create a new configuration
    else
        call iniconfig(x, y, z, d)
    end if
    ! Energy of the initial configuration
    call force <<< grid, tBlock >>> (x, y, z, fx_d, fy_d, fz_d, enerpot, zfac)
    istat = cudaDeviceSynchronize()
    if ( istat .ne. 0 ) then
        write(unit=u, fmt='(a)') 'Error with GPU!'
    end if
    write(unit=output_unit, fmt='(a,8f)') 'E/N=', sum(enerpot) / real(np)

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

    write(unit=output_unit, fmt='(a)') 'The system has thermalized'

    ! Save positions to file for re-use in other runs
    open(newunit=u, file='finalconBD.dat', status='unknown')
    do i = 1, np
        write(u, '(3f16.8)') x(i), y(i), z(i)
    end do
    close(u)

    ! Open file for production run
    open(newunit=u, file='energy_BD_prod.dat', status='unknown')
    
    ! Turn off periodic boundary conditions
    pbc = 0
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
            ! Save snapshots to file
            call snapshots(x, y, z, istep, 'production.xyz')
        end if
    end do
    
    ! Close file for production
    close(u)
end program main
