module utils
    use parameters
    
    implicit none
    
    public iniconfig

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
end module utils