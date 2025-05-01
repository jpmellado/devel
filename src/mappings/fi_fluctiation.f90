!########################################################################
! Reynolds fluctuations of array a
!########################################################################
subroutine FI_FLUCTUATION_INPLACE(nx, ny, nz, a)
    use TLab_Constants, only: wp, wi
    use Averages, only: AVG_IK
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: a(nx, ny, nz)

    ! -------------------------------------------------------------------
    real(wp) dummy
    integer(wi) k

    ! ###################################################################
    do k = 1, nz
        dummy = AVG_IK(nx, ny, nz, k, a)
        a(:, :, k) = a(:, :, k) - dummy
    end do

    return
end subroutine FI_FLUCTUATION_INPLACE
