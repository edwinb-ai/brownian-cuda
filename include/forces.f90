module forcec
    use parameters
    use cudafor

    implicit none

    ! Intrinsic pseudohs potential arguments that need to be computed
    ! only once
    real, parameter :: dlr = 50.0
    real, parameter :: dT = 1.4737
    real, parameter :: dla = 49.0
    real, parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
    real, parameter :: rchs = (dlr/dla)**(1.0/(dlr-dla))

    public force
    
contains

    attributes(global) subroutine force(x, y, z, fx, fy, fz, ener, zfac)
        real, intent(in) :: x(:), y(:), z(:)
        real, intent(inout) :: fx(:), fy(:), fz(:)
        real, intent(inout) :: ener(:), zfac(:)

        ! Local variables
        integer :: i, j, s, idx
        real :: xij, yij, zij, uij, fij
        real :: potential
        real :: fxij, fyij, fzij, rij
        real :: fxx, fyy, fzz, virial, vtmp, rinv
    
        ! Inicialización de variables de CUDA
        idx = blockDim%x * (blockIdx%x - 1) + threadIdx%x
        s = blockDim%x * gridDim%x
        
        do i = idx, np, s
            potential = 0.0
            virial = 0.0
            fxx = 0.0
            fyy = 0.0
            fzz = 0.0

            do j = 1, np
                if ( i == j ) cycle

                xij = x(i) - x(j)
                yij = y(i) - y(j)
                zij = z(i) - z(j)

                xij = xij - boxl*nint(xij / boxl)
                yij = yij - boxl*nint(yij / boxl)
                zij = zij - boxl*nint(zij / boxl)
                
                rij = sqrt(xij*xij + yij*yij + zij*zij)
                
                if ( (rij > 0.1) .and. (rij < rc) ) then
                    rinv = 1.0 / rij
                    if ( rij < rchs ) then
                        uij = a2 * ((rinv**dlr)-(rinv**dla))
                        uij = uij + 1.0
                        fij = (dlr*(rinv**(dlr+1.0))) - (dla * (rinv**dlr))
                        fij = fij * a2
                    else
                        uij = 0.0
                        fij = 0.0
                    end if

                    fxij = fij * xij * rinv
                    fyij = fij * yij * rinv
                    fzij = fij * zij * rinv

                    ! Accumulate energy and virial
                    potential = potential + uij
                    vtmp = (fxij * xij) + (fyij * yij) + (fzij * zij)
                    virial = virial + vtmp

                    fxx = fxx + fxij
                    fyy = fyy + fyij
                    fzz = fzz + fzij
                end if
            end do
            
            fx(i) = fxx
            fy(i) = fyy
            fz(i) = fzz

            ener(i) = potential / 2.0
            zfac(i) = virial / 2.0
        end do
    end subroutine force
end module