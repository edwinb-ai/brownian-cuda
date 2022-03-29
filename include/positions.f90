module positions
    use parameters
    use cudafor

    implicit none
    
    public position

    contains

    attributes(global) subroutine position(x,y,z,fx,fy,fz,randv,pbc)
        real, intent(in) :: fx(:), fy(:), fz(:)
        real, intent(inout) :: x(:), y(:), z(:)
        real, intent(inout) :: randv(:)
        integer, intent(in), value :: pbc
        ! Local variables
        integer :: i, s, idx

        idx = blockDim%x * (blockIdx%x - 1) + threadIdx%x
        s = blockDim%x * gridDim%x

        do i = idx, np, s
            x(i) = x(i) + (fx(i) * deltat / ktemp) + randv(i)
            y(i) = y(i) + (fy(i) * deltat / ktemp) + randv(i+np)
            z(i) = z(i) + (fz(i) * deltat / ktemp) + randv(i+(2*np))

            if ( pbc == 1 ) then
                x(i) = x(i) - boxl*nint(x(i)/boxl)
                y(i) = y(i) - boxl*nint(y(i)/boxl)
                z(i) = z(i) - boxl*nint(z(i)/boxl)
            end if
        end do
    end subroutine position
end module