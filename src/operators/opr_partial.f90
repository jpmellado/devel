module OPR_Partial
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_NONE
    use TLab_Arrays, only: wrk2d, wrk3d
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
#endif
    use FDM, only: fdm_dt
    use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
    implicit none
    private

    public :: OPR_Partial_Initialize
    public :: OPR_Partial_X     ! These first 2 could be written in terms of the last one...
    public :: OPR_Partial_Y
    public :: OPR_Partial_Z

    integer, parameter, public :: OPR_P1 = 1                ! 1. order derivative
    integer, parameter, public :: OPR_P2 = 2                ! 2. order derivative
    integer, parameter, public :: OPR_P2_P1 = 3             ! 2. and 1.order derivatives

    ! -----------------------------------------------------------------------
    procedure(OPR_Partial_interface) :: OPR_Partial_dt
    abstract interface
        subroutine OPR_Partial_interface(type, nx, ny, nz, g, u, result, tmp1, ibc)
            use TLab_Constants, only: wi, wp
            use FDM, only: fdm_dt
            integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
            integer(wi), intent(in) :: nx, ny, nz
            type(fdm_dt), intent(in) :: g
            real(wp), intent(in) :: u(nx*ny*nz)
            real(wp), intent(out) :: result(nx*ny*nz)
            real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
            integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative
        end subroutine
    end interface
    procedure(OPR_Partial_dt), pointer :: OPR_Partial_X, OPR_Partial_Y

contains
    ! ###################################################################
    ! ###################################################################
    subroutine OPR_Partial_Initialize()

        ! ###################################################################
        ! Setting procedure pointers
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            OPR_Partial_X => OPR_Partial_X_Parallel
        else
#endif
            OPR_Partial_X => OPR_Partial_X_Serial
#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (ims_npro_j > 1) then
            OPR_Partial_Y => OPR_Partial_Y_Parallel
        else
#endif
            OPR_Partial_Y => OPR_Partial_Y_Serial
#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Partial_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine OPR_Partial_X_Serial(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer ibc_loc

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, g%size, g%size, ny*nz, result, ny*nz)
#else
        call TLab_Transpose(u, g%size, ny*nz, g%size, result, ny*nz)
#endif

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(ny*nz, ibc_loc, g%der1, g%der1%lu, result, tmp1, wrk2d)
            call FDM_Der2_Solve(ny*nz, g%der2, g%der2%lu, result, wrk3d, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(ny*nz, ibc_loc, g%der1, g%der1%lu, result, wrk3d, wrk2d)
            call FDM_Der2_Solve(ny*nz, g%der2, g%der2%lu, result, tmp1, wrk3d, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(ny*nz, ibc_loc, g%der1, g%der1%lu, result, wrk3d, wrk2d)

        end select

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, ny*nz, ny*nz, g%size, result, g%size)
            call DGETMO(wrk3d, ny*nz, ny*nz, g%size, tmp1, g%size)
        else
            call DGETMO(wrk3d, ny*nz, ny*nz, g%size, result, g%size)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, ny*nz, g%size, ny*nz, result, g%size)
            call TLab_Transpose(wrk3d, ny*nz, g%size, ny*nz, tmp1, g%size)
        else
            call TLab_Transpose(wrk3d, ny*nz, g%size, ny*nz, result, g%size)
        end if
#endif

        return
    end subroutine OPR_Partial_X_Serial

    ! ###################################################################
    ! ###################################################################
#ifdef USE_MPI
    subroutine OPR_Partial_X_Parallel(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer ibc_loc

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        nlines = tmpi_plan_dx%nlines

        ! Transposition: make x-direction the last one
        call TLabMPI_Trp_ExecI_Forward(u, result, tmpi_plan_dx)
#ifdef USE_ESSL
        call DGETMO(result, g%size, g%size, nlines, wrk3d, nlines)
#else
        call TLab_Transpose(result, g%size, nlines, g%size, wrk3d, nlines)
#endif

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, wrk3d, result, wrk2d)

        end select

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, nlines, nlines, g%size, wrk3d, g%size)
            call DGETMO(result, nlines, nlines, g%size, tmp1, g%size)
            call TLabMPI_Trp_ExecI_Backward(tmp1, result, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, tmp1, tmpi_plan_dx)
        else
            call DGETMO(result, nlines, nlines, g%size, wrk3d, g%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, nlines, g%size, nlines, wrk3d, g%size)
            call TLab_Transpose(result, nlines, g%size, nlines, tmp1, g%size)
            call TLabMPI_Trp_ExecI_Backward(tmp1, result, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, tmp1, tmpi_plan_dx)
        else
            call TLab_Transpose(result, nlines, g%size, nlines, wrk3d, g%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
        end if
#endif

        return
    end subroutine OPR_Partial_X_Parallel
#endif

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_Y_Serial(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer ibc_loc

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose(u, nx*ny, nz, nx*ny, result, nz)
#endif
        nlines = nx*nz

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, result, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, result, wrk3d, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, result, wrk3d, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, result, tmp1, wrk3d, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, result, wrk3d, wrk2d)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, nz, nz, nx*ny, result, nx*ny)
            call DGETMO(wrk3d, nz, nz, nx*ny, tmp1, nx*ny)
        else
            call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, nz, nx*ny, nz, result, nx*ny)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)
        else
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
        end if
#endif

        return
    end subroutine OPR_Partial_Y_Serial

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine OPR_Partial_Y_Parallel(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer ibc_loc

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose(u, nx*ny, nz, nx*ny, result, nz)
#endif
        call TLabMPI_Trp_ExecJ_Forward(result, wrk3d, tmpi_plan_dy)
        nlines = tmpi_plan_dy%nlines

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g%der2, g%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g%der1, g%der1%lu, wrk3d, result, wrk2d)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call TLabMPI_Trp_ExecJ_Backward(tmp1, wrk3d, tmpi_plan_dy)
            call TLabMPI_Trp_ExecJ_Backward(result, tmp1, tmpi_plan_dy)
            call DGETMO(tmp1, nz, nz, nx*ny, result, nx*ny)
            call DGETMO(wrk3d, nz, nz, nx*ny, tmp1, nx*ny)
        else
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLabMPI_Trp_ExecJ_Backward(tmp1, wrk3d, tmpi_plan_dy)
            call TLabMPI_Trp_ExecJ_Backward(result, tmp1, tmpi_plan_dy)
            call TLab_Transpose(tmp1, nz, nx*ny, nz, result, nx*ny)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)
        else
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
        end if
#endif

        return
    end subroutine OPR_Partial_Y_Parallel
#endif

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_Z(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer ibc_loc

        ! ###################################################################
        if (g%size == 1) then
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(nx*ny, ibc_loc, g%der1, g%der1%lu, u, wrk3d, wrk2d)
            call FDM_Der2_Solve(nx*ny, g%der2, g%der2%lu, u, result, wrk3d, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nx*ny, ibc_loc, g%der1, g%der1%lu, u, tmp1, wrk2d)
            call FDM_Der2_Solve(nx*ny, g%der2, g%der2%lu, u, result, tmp1, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nx*ny, ibc_loc, g%der1, g%der1%lu, u, result, wrk2d)

        end select

        return
    end subroutine OPR_Partial_Z

end module OPR_Partial
