module utils
    use parameters
    
    implicit none
    
    public iniconfig, snapshots, parse_input
contains
    subroutine iniconfig(xc, yc, zc, d)
    ! defining three vector of mp dimension, it indicate that only are out variables
        real, intent(out) :: xc(:), yc(:), zc(:)
        real, intent(in) :: d
        ! Local variables
        integer :: i

        xc(1) = -1.0*(boxl-d) / 2.0
        yc(1) = -1.0*(boxl-d) / 2.0
        zc(1) = -1.0*(boxl-d) / 2.0

        do i = 2,np
            xc(i) = xc(i-1) + d
            yc(i) = yc(i-1)
            zc(i) = zc(i-1)
            if (xc(i) > rc) then
                xc(i) = xc(1)
                yc(i) = yc(i-1) + d
                if (yc(i) > rc) then
                    xc(i) = xc(1)
                    yc(i) = yc(1)
                    zc(i) = zc(i-1) + d
                end if
            end if
        end do
    end subroutine iniconfig

    subroutine parse_input(filein, limg, limt, initconf)
        character(len=*), intent(in) :: filein
        logical, intent(inout) :: initconf
        integer, intent(inout) :: limg, limt

        !! Local variables
        integer :: u

        open(newunit=u, file=filein, status="old", action="read")
        !! Read the inputs from the file
        read(u, *) phi ! Density
        read(u, *) limt ! Thermalization steps
        read(u, *) limg ! Averaging steps
        read(u, *) initconf
        close(u)
        
    end subroutine parse_input

    subroutine snapshots(x, y, z, istep, filename)
        ! Arguments
        real, intent(in) :: x(:), y(:), z(:)
        integer, intent(in) :: istep
        character(len=*), intent(in) :: filename

        ! Local variables
        logical :: exists
        integer :: io, i, arrsize

        arrsize = size(x)

        inquire(file=filename, exist=exists)
        if (exists) then
            open(newunit=io, file=filename, position="append", &
                & status="old", action="write")
            write(io, fmt='(i4)') arrsize
            write(io, fmt='(i6.4)') istep
        else
            open(newunit=io, file=filename, status="new", action="write")
            write(io, fmt='(i4)') arrsize
            write(io, fmt='(i6.4)') istep
        end if

        do i = 1, arrsize
            write(io, fmt='(i1,A,f11.8,A,f12.8,A,f12.8)') 1, ' ', x(i), ' ', &
                & y(i), ' ', z(i)
        end do

        close(io)
    end subroutine snapshots
end module utils