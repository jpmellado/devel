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
    ! ! #######################################################################
    ! ! #######################################################################
    ! subroutine Reduce_Block_InPlace(nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst, a, wrk1d)
    !     use TLab_Constants, only: wp, wi
    !     implicit none

    !     integer(wi), intent(IN) :: nx, ny, nz, nx1, ny1, nz1, nx_dst, ny_dst, nz_dst
    !     real(wp), dimension(nx*ny*nz), intent(INOUT) :: a
    !     real(wp), dimension(nx_dst), intent(INOUT) :: wrk1d

    ! ! -------------------------------------------------------------------
    !     integer(wi) j, k, nxy, nxy_dst, ip, ip_dst

    ! ! -------------------------------------------------------------------
    !     nxy = nx*ny
    !     nxy_dst = nx_dst*ny_dst

    !     do k = 1, nz_dst
    !         ip = (k - 1)*nxy + (nz1 - 1)*nxy + (ny1 - 1)*nx + (nx1 - 1) + 1
    !         ip_dst = (k - 1)*nxy_dst + 1
    !         do j = 1, ny_dst
    !             wrk1d(1:nx_dst) = a(ip:ip + nx_dst - 1)
    !             a(ip_dst:ip_dst + nx_dst - 1) = wrk1d(1:nx_dst)

    !             ip = ip + nx
    !             ip_dst = ip_dst + nx_dst
    !         end do
    !     end do

    !     return
    ! end subroutine Reduce_Block_InPlace

end module Reductions
