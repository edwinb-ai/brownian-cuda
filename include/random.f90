module randomf
    use curand
    use parameters
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    real, parameter, private :: sqrt3 = sqrt(3.0)
    
    public initialize_rng, fill_vec
contains
    subroutine initialize_rng(rng)
        type(curandGenerator), intent(inout) :: rng

        ! Local variables
        real :: some_seed
        integer, allocatable :: seed(:)
        integer :: host_seed, istat, largest, io, n

        ! Initialize PRNG
        call random_seed(size=n)
        allocate(seed(n))
        ! Read from urandom file
        open(newunit=io, file="/dev/urandom", access="stream", form="unformatted", &
            action="read", status="old", iostat=istat)
        if (istat==0) then
            read(io, iostat=istat) seed
            close(io)
        endif
        ! Set the final seed
        call random_seed(put=seed)
        ! Initialize the PRNG for the device
        istat = curandCreateGenerator(rng, CURAND_RNG_PSEUDO_DEFAULT)
        ! Compute the seed for the largest 32-bit integer
        call random_number(some_seed)
        host_seed = floor(some_seed * huge(largest))
        ! Show the seed to screen
        write(unit=output_unit, fmt='(a,i)') 'Random seed =', host_seed
        istat = curandSetPseudoRandomGeneratorSeed(rng, host_seed)
        if ( istat .ne. 0 ) then
            write(unit=output_unit, fmt='(a)') 'Error with GPU!'
        end if

    end subroutine initialize_rng

    subroutine fill_vec(vec, rng, threenp)
        real, device, intent(inout) :: vec(:)
        type(curandGenerator), intent(in) :: rng
        integer, intent(in) :: threenp
    
        ! Local variables
        integer :: istat, i

        ! Initialize the vectors using a uniform distribution
        istat = curandGenerateUniform(rng, vec, threenp)
        if ( istat .ne. 0 ) then
            write(unit=output_unit, fmt='(a)') 'Error with GPU!'
        end if
        
        ! Move the distribution to (-1, 1)
        ! Multiply by sqrt(6.0 * deltat) to make the second moment
        ! equal as in a normal distribution

        !$cuf kernel do <<< *,* >>>
        do i = 1, threenp
            vec(i) = sqrt3 * sqtwodt * ((2.0 * vec(i)) - 1.0)
        end do
    end subroutine fill_vec

    subroutine fill_vec_normal(vec, rng, threenp)
        real, device, intent(inout) :: vec(:)
        type(curandGenerator), intent(in) :: rng
        integer, intent(in) :: threenp
    
        ! Local variables
        integer :: istat, i

        ! Initialize the vectors using a uniform distribution
        istat = curandGenerateNormal(rng, vec, threenp, 0.0, sqtwodt)
        if ( istat .ne. 0 ) then
            write(unit=output_unit, fmt='(a)') 'Error with GPU!'
        end if
    end subroutine fill_vec_normal
end module randomf