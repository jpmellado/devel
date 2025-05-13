module Reductions
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: Reduce_Block_InPlace

contains
    ! #######################################################################
    ! #######################################################################
    subroutine Reduce_Block_InPlace(nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst, a)
        integer(wi), intent(in) :: nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst
        real(wp), intent(inout), target :: a(nx, ny, nz)

        integer(wi) i, j, k
        real(wp), pointer :: a_dst(:, :, :) => null()

        a_dst(1:nx_dst, 1:ny_dst, 1:nz_dst) => a(1:nx_dst*ny_dst*nz_dst, 1, 1)

        do k = 1, nz_dst
            do j = 1, ny_dst
                do i = 1, nx_dst
                    a_dst(i, j, k) = a(nx1 + i - 1, ny1 + j - 1, nz1 + k - 1)
                end do
            end do
        end do

        return
    end subroutine Reduce_Block_InPlace

end module Reductions
