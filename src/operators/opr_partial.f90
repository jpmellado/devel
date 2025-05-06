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
    ! use FDM_Interpolate, only: FDM_Interpol, FDM_Interpol_Der1
    implicit none
    private

    public :: OPR_Partial_X     ! These first 2 could be written in terms of the last one...
    public :: OPR_Partial_Y
    public :: OPR_Partial_Z

    integer, parameter, public :: OPR_P1 = 1                ! 1. order derivative
    integer, parameter, public :: OPR_P2 = 2                ! 2. order derivative
    integer, parameter, public :: OPR_P2_P1 = 3             ! 2. and 1.order derivatives

contains
    ! ###################################################################
    ! ###################################################################
    subroutine OPR_Partial_X(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        target u, tmp1, result

        ! -------------------------------------------------------------------
        integer(wi) nyz
        integer ibc_loc

        real(wp), dimension(:), pointer :: p_a, p_b, p_c, p_d

        ! ###################################################################
        if (g%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call TLabMPI_Trp_ExecI_Forward(u, result, tmpi_plan_dx)
            p_a => result
            p_b => wrk3d
            p_c => result
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_d => tmp1
            end if
            nyz = tmpi_plan_dx%nlines
        else
#endif
            p_a => u
            p_b => result
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_c => tmp1
                p_d => wrk3d
            else
                p_c => wrk3d
            end if
            nyz = ny*nz
#ifdef USE_MPI
        end if
#endif

#ifdef USE_ESSL
        call DGETMO(p_a, g%size, g%size, nyz, p_b, nyz)
#else
        call TLab_Transpose(p_a, g%size, nyz, g%size, p_b, nyz)
#endif

        ! ###################################################################
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(nyz, ibc_loc, g%der1, g%der1%lu, p_b, p_d, wrk2d)
            call FDM_Der2_Solve(nyz, g%der2, g%der2%lu, p_b, p_c, p_d, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nyz, ibc_loc, g%der1, g%der1%lu, p_b, p_d, wrk2d)
            call FDM_Der2_Solve(nyz, g%der2, g%der2%lu, p_b, p_c, p_d, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nyz, ibc_loc, g%der1, g%der1%lu, p_b, p_c, wrk2d)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(p_c, nyz, nyz, g%size, p_b, g%size)
        if (type == OPR_P2_P1) call DGETMO(p_d, nyz, nyz, g%size, p_c, g%size)
#else
        call TLab_Transpose(p_c, nyz, g%size, nyz, p_b, g%size)
        if (type == OPR_P2_P1) call TLab_Transpose(p_d, nyz, g%size, nyz, p_c, g%size)
#endif

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            if (type == OPR_P2_P1) call TLabMPI_Trp_ExecI_Backward(p_c, tmp1, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(p_b, result, tmpi_plan_dx)
        end if
#endif

        nullify (p_a, p_b, p_c)
        if (associated(p_d)) nullify (p_d)

        return
    end subroutine OPR_Partial_X

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_Y(type, nx, ny, nz, g, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        target u, tmp1, result

        ! -------------------------------------------------------------------
        integer(wi) nxz
        integer ibc_loc
        real(wp), dimension(:), pointer :: p_b, p_c, p_d

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

#ifdef USE_MPI
        if (ims_npro_j > 1) then
            call TLabMPI_Trp_ExecJ_Forward(result, wrk3d, tmpi_plan_dy)
            nxz = tmpi_plan_dy%nlines

            p_b => wrk3d
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_c => result
                p_d => tmp1
            else
                p_c => result
            end if

        else
#endif
            nxz = nx*nz

            p_b => result
            if (any([OPR_P2, OPR_P2_P1] == type)) then
                p_c => tmp1
                p_d => wrk3d
            else
                p_c => wrk3d
            end if

#ifdef USE_MPI
        end if
#endif

        ! ###################################################################
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g%der2%need_1der) call FDM_Der1_Solve(nxz, ibc_loc, g%der1, g%der1%lu, p_b, p_d, wrk2d)
            call FDM_Der2_Solve(nxz, g%der2, g%der2%lu, p_b, p_c, p_d, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nxz, ibc_loc, g%der1, g%der1%lu, p_b, p_d, wrk2d)
            call FDM_Der2_Solve(nxz, g%der2, g%der2%lu, p_b, p_c, p_d, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nxz, ibc_loc, g%der1, g%der1%lu, p_b, p_c, wrk2d)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_MPI
        if (ims_npro_j > 1) then
            call TLabMPI_Trp_ExecJ_Backward(p_c, p_b, tmpi_plan_dy)
            if (type == OPR_P2_P1) call TLabMPI_Trp_ExecJ_Backward(p_d, p_c, tmpi_plan_dy)

#ifdef USE_ESSL
            if (type == OPR_P2_P1) call DGETMO(p_c, nz, nz, nx*ny, tmp1, nx*ny)
            call DGETMO(p_b, nz, nz, nx*ny, result, nx*ny)
#else
            if (type == OPR_P2_P1) call TLab_Transpose(p_c, nz, nx*ny, nz, tmp1, nx*ny)
            call TLab_Transpose(p_b, nz, nx*ny, nz, result, nx*ny)
#endif

        else
#endif

#ifdef USE_ESSL
            call DGETMO(p_c, nz, nz, nx*ny, result, nx*ny)
            if (type == OPR_P2_P1) call DGETMO(p_d, nz, nz, nx*ny, tmp1, nx*ny)
#else
            call TLab_Transpose(p_c, nz, nx*ny, nz, result, nx*ny)
            if (type == OPR_P2_P1) call TLab_Transpose(p_d, nz, nx*ny, nz, tmp1, nx*ny)
#endif

#ifdef USE_MPI
        end if
#endif

        nullify (p_b, p_c)
        if (associated(p_d)) nullify (p_d)

        return
    end subroutine OPR_Partial_Y

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
