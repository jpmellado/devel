#include "tlab_error.h"

! Calculate the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s

module NSE_Burgers
    use TLab_Constants, only: wp, wi, efile, lfile, BCS_NONE
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk2d, wrk3d
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: fdm_dt, g
    use NavierStokes, only: visc, schmidt
    use OPR_Partial
    use LargeScaleForcing, only: subsidenceProps, TYPE_SUB_CONSTANT, wbackground
    implicit none
    private

    public :: NSE_Burgers_Initialize
    public :: NSE_Burgers_X
    public :: NSE_Burgers_Y
    public :: NSE_Burgers_Z

    ! -----------------------------------------------------------------------
    procedure(NSE_Burgers_interface) :: NSE_Burgers_dt
    abstract interface
        subroutine NSE_Burgers_interface(is, nx, ny, nz, s, result, tmp1, u_t)
            use TLab_Constants, only: wi, wp
            integer, intent(in) :: is                       ! scalar index; if 0, then velocity
            integer(wi), intent(in) :: nx, ny, nz
            real(wp), intent(in) :: s(nx*ny*nz)
            real(wp), intent(out) :: result(nx*ny*nz)
            real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
            real(wp), intent(in), optional, target :: u_t(nx*ny*nz)
        end subroutine
    end interface
    procedure(NSE_Burgers_dt), pointer :: NSE_Burgers_X, NSE_Burgers_Y

    type :: fdm_diffusion_dt
        sequence
        real(wp), allocatable :: lu(:, :, :)
        ! type(rho_anelastic_dt) :: rho_anelastic
    end type fdm_diffusion_dt
    type(fdm_diffusion_dt) :: fdmDiffusion(3)

    type :: rho_anelastic_dt                        ! 1/rho in diffusion term in anelastic formulation
        sequence
        logical :: active = .false.
        real(wp), allocatable :: values(:)
    end type rho_anelastic_dt
    type(rho_anelastic_dt) :: rho_anelastic(3)      ! one for each direction

    real(wp), dimension(:), pointer :: p_vel

    ! type(filter_dt) :: Dealiasing(3)
    integer :: Dealiasing(3) ! tobefixed
    ! real(wp), allocatable, target :: wrkdea(:, :)       ! Work arrays for dealiasing (scratch space)

contains
    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_Memory, only: inb_scal
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_j, ims_npro_j
#endif
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
        use Thermo_Anelastic, only: rbackground, ribackground

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile

        integer(wi) ig, is, ip, j, idummy
        integer(wi) nlines, offset
        real(wp) dummy

        ! ###################################################################
        ! Read input data
        bakfile = trim(adjustl(inifile))//'.bak'

        ! call FILTER_READBLOCK(bakfile, inifile, 'Dealiasing', Dealiasing)

        ! ###################################################################
        ! Initialize LU factorization of the second-order derivative times the diffusivity
        do ig = 1, 3
            if (g(ig)%size == 1) cycle

            if (g(ig)%der2%nb_diag(1) /= 3) then
                call TLab_Write_ASCII(efile, __FILE__//'. Undeveloped for more than 3 LHS diagonals in 2. order derivatives.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            if (g(ig)%periodic) then
                idummy = g(ig)%der2%nb_diag(1) + 2
            else
                idummy = g(ig)%der2%nb_diag(1)
            end if
            allocate (fdmDiffusion(ig)%lu(g(ig)%size, idummy, 0:inb_scal))

            ip = 0
            do is = 0, inb_scal ! case 0 for the reynolds number
                if (is == 0) then
                    dummy = visc
                else
                    dummy = visc/schmidt(is)
                end if

                if (g(ig)%periodic) then                        ! Check routines TRIDPFS and TRIDPSS
                    fdmDiffusion(ig)%lu(:, 1, is) = g(ig)%der2%lu(:, 1)         ! matrix L; 1. subdiagonal
                    fdmDiffusion(ig)%lu(:, 2, is) = g(ig)%der2%lu(:, 2)*dummy   ! matrix L; 1/diagonal
                    fdmDiffusion(ig)%lu(:, 3, is) = g(ig)%der2%lu(:, 3)         ! matrix U is the same
                    fdmDiffusion(ig)%lu(:, 4, is) = g(ig)%der2%lu(:, 4)/dummy   ! matrix L; Additional row/column
                    fdmDiffusion(ig)%lu(:, 5, is) = g(ig)%der2%lu(:, 5)         ! matrix U is the same

                else                                            ! Check routines TRIDFS and TRIDSS
                    fdmDiffusion(ig)%lu(:, 1, is) = g(ig)%der2%lu(:, 1)         ! matrix L is the same
                    fdmDiffusion(ig)%lu(:, 2, is) = g(ig)%der2%lu(:, 2)*dummy   ! matrix U; 1/diagonal
                    fdmDiffusion(ig)%lu(:, 3, is) = g(ig)%der2%lu(:, 3)/dummy   ! matrix U; 1. superdiagonal

                end if

            end do
        end do

        ! ###################################################################
        ! Initialize dealiasing
        ! do ig = 1, 3
        !     if (Dealiasing(ig)%type /= DNS_FILTER_NONE) call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
        ! end do

        ! if (any(Dealiasing(:)%type /= DNS_FILTER_NONE)) then
        !     call TLab_Allocate_Real(__FILE__, wrkdea, [isize_field, 2], 'wrk-dealiasing')
        ! end if

        ! ###################################################################
        ! Initialize anelastic density correction
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call TLab_Write_ASCII(lfile, 'Initialize anelastic density correction in burgers operator.')

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along X
            rho_anelastic(1)%active = .true.
#ifdef USE_MPI
            if (ims_npro_i > 1) then
                nlines = tmpi_plan_dx%nlines
                offset = nlines*ims_pro_i
            else
#endif
                nlines = jmax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho_anelastic(1)%values(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, y%size) + 1
                rho_anelastic(1)%values(j) = ribackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Y
            rho_anelastic(2)%active = .true.
#ifdef USE_MPI
            if (ims_npro_j > 1) then
                nlines = tmpi_plan_dy%nlines
                offset = nlines*ims_pro_j
            else
#endif
                nlines = imax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho_anelastic(2)%values(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, z%size) + 1
                rho_anelastic(2)%values(j) = ribackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Z; see FDM_CreatePlan
            ! we implement it directly in the tridiagonal system
            do is = 0, inb_scal ! case 0 for the velocity
                fdmDiffusion(3)%lu(:, 2, is) = fdmDiffusion(3)%lu(:, 2, is)*ribackground(:)  ! matrix U; 1/diagonal
                fdmDiffusion(3)%lu(:z%size - 1, 3, is) = fdmDiffusion(3)%lu(:z%size - 1, 3, is)*rbackground(2:) ! matrix U; 1. superdiagonal
            end do

        end if

        ! ###################################################################
        ! Setting procedure pointers
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            NSE_Burgers_X => NSE_Burgers_X_Parallel
        else
#endif
            NSE_Burgers_X => NSE_Burgers_X_Serial
#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (ims_npro_j > 1) then
            NSE_Burgers_Y => NSE_Burgers_Y_Parallel
        else
#endif
            NSE_Burgers_Y => NSE_Burgers_Y_Serial
#ifdef USE_MPI
        end if
#endif
        return
    end subroutine NSE_Burgers_Initialize

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_X_Serial(is, nx, ny, nz, s, result, tmp1, u_t)
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
        real(wp), intent(in), optional, target :: u_t(nx*ny*nz)

        ! -------------------------------------------------------------------
        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, g(1)%size, g(1)%size, ny*nz, tmp1, ny*nz)
#else
        call TLab_Transpose(s, g(1)%size, ny*nz, g(1)%size, tmp1, ny*nz)
#endif

        if (present(u_t)) then  ! transposed velocity is passed as argument
            p_vel => u_t
        else
            p_vel => tmp1
        end if

        call NSE_Burgers_1D(ny*nz, g(1), fdmDiffusion(1)%lu(:, :, is), rho_anelastic(1), Dealiasing(1), tmp1, p_vel, wrk3d, result)

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, g(1)%size, result, g(1)%size)
#else
        call TLab_Transpose(wrk3d, ny*nz, g(1)%size, ny*nz, result, g(1)%size)
#endif

        return
    end subroutine NSE_Burgers_X_Serial

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_Burgers_X_Parallel(is, nx, ny, nz, s, result, tmp1, u_t)
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
        real(wp), intent(in), optional, target :: u_t(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        nlines = tmpi_plan_dx%nlines

        ! Transposition: make x-direction the last one
        call TLabMPI_Trp_ExecI_Forward(s, result, tmpi_plan_dx)
#ifdef USE_ESSL
        call DGETMO(result, g(1)%size, g(1)%size, nlines, tmp1, nlines)
#else
        call TLab_Transpose(result, g(1)%size, nlines, g(1)%size, tmp1, nlines)
#endif

        if (present(u_t)) then  ! transposed velocity is passed as argument
            p_vel => u_t
        else
            p_vel => tmp1
        end if

        call NSE_Burgers_1D(nlines, g(1), fdmDiffusion(1)%lu(:, :, is), rho_anelastic(1), Dealiasing(1), tmp1, p_vel, result, wrk3d)

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(result, nlines, nlines, g(1)%size, wrk3d, g(1)%size)
#else
        call TLab_Transpose(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
#endif
        call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)

        return
    end subroutine NSE_Burgers_X_Parallel
#endif

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_Y_Serial(is, nx, ny, nz, s, result, tmp1, u_t)
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
        real(wp), intent(in), target, optional :: u_t(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx*ny, nx*ny, nz, tmp1, nz)
#else
        call TLab_Transpose(s, nx*ny, nz, nx*ny, tmp1, nz)
#endif
        nlines = nx*nz

        if (present(u_t)) then  ! transposed velocity is passed as argument
            p_vel => u_t
        else
            p_vel => tmp1
        end if

        call NSE_Burgers_1D(nlines, g(2), fdmDiffusion(2)%lu(:, :, is), rho_anelastic(2), Dealiasing(2), tmp1, p_vel, wrk3d, result)

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
#else
        call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
#endif

        return
    end subroutine NSE_Burgers_Y_Serial

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_Burgers_Y_Parallel(is, nx, ny, nz, s, result, tmp1, u_t)
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
        real(wp), intent(in), target, optional :: u_t(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx*ny, nx*ny, nz, wrk3d, nz)
#else
        call TLab_Transpose(s, nx*ny, nz, nx*ny, wrk3d, nz)
#endif
        call TLabMPI_Trp_ExecJ_Forward(wrk3d, tmp1, tmpi_plan_dy)
        nlines = tmpi_plan_dy%nlines

        if (present(u_t)) then  ! transposed velocity is passed as argument
            p_vel => u_t
        else
            p_vel => tmp1
        end if

        call NSE_Burgers_1D(nlines, g(2), fdmDiffusion(2)%lu(:, :, is), rho_anelastic(2), Dealiasing(2), tmp1, p_vel, result, wrk3d)

        ! Put arrays back in the order in which they came in
        call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
#else
        call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
#endif

        return
    end subroutine NSE_Burgers_Y_Parallel
#endif

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_Z(is, nx, ny, nz, s, result, u)
        use TLab_Pointers_2D, only: p2d_wrk3d
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny, nz)
        real(wp), intent(in) :: u(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) k

        ! ###################################################################
        if (z%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        call NSE_Burgers_1D(nx*ny, g(3), fdmDiffusion(3)%lu(:, :, is), rho_anelastic(3), Dealiasing(3), s, u, result, wrk3d)

        if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
            do k = 1, nz
                result(:, k) = result(:, k) + p2d_wrk3d(:, k)*wbackground(k)
            end do
        end if

        return
    end subroutine NSE_Burgers_Z

    !########################################################################
    !# Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s
    !# along generic direction x to nlines lines of data
    !#
    !# Second derivative uses LE decomposition including diffusivity coefficient
    !########################################################################
    subroutine NSE_Burgers_1D(nlines, g, lu2d, rhoi, dealiasing, s, u, result, dsdx)
        use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
        integer(wi), intent(in) :: nlines       ! # of lines to be solved
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: lu2d(:, :)      ! LU decomposition including the diffusion parameter for corresponding field is
        type(rho_anelastic_dt), intent(in) :: rhoi
        ! type(filter_dt), intent(in) :: dealiasing
        integer :: dealiasing
        real(wp), intent(in) :: s(nlines, g%size), u(nlines, g%size)  ! argument field and velocity field
        real(wp), intent(out) :: result(nlines, g%size)                ! N(u) applied to s
        real(wp), intent(inout) :: dsdx(nlines, g%size)                  ! dsdx

        ! -------------------------------------------------------------------
        integer(wi) ij
        ! real(wp), pointer :: uf(:, :) !, dsf(:, :)

        ! ###################################################################
        ! dsdx: 1st derivative; result: 2nd derivative including diffusivity
        call FDM_Der1_Solve(nlines, BCS_NONE, g%der1, g%der1%lu, s, dsdx, wrk2d)
        call FDM_Der2_Solve(nlines, g%der2, lu2d, s, result, dsdx, wrk2d)

        ! ###################################################################
        ! Operation; diffusivity included in 2.-order derivativelu2_p
        ! ###################################################################
        ! if (dealiasing%type /= DNS_FILTER_NONE) then
        !     uf(1:nlines, 1:g%size) => wrkdea(1:nlines*g%size, 1)
        !     dsf(1:nlines, 1:g%size) => wrkdea(1:nlines*g%size, 2)
        !     call OPR_FILTER_1D(nlines, dealiasing, u, uf)
        !     call OPR_FILTER_1D(nlines, dealiasing, dsdx, dsf)

        !     result(:, :) = result(:, :) - uf(:, :)*dsf(:, :)

        !     nullify (uf, dsf)

        ! else
        ! result(:, :) = result(:, :) - u(:, :)*dsdx(:, :)
        ! end if
        if (rhoi%active) then
            do ij = 1, g%size
                result(:, ij) = result(:, ij)*rhoi%values(:) - u(:, ij)*dsdx(:, ij)
            end do

        else
            result(:, :) = result(:, :) - u(:, :)*dsdx(:, :)
        end if

        return
    end subroutine NSE_Burgers_1D

end module NSE_Burgers
