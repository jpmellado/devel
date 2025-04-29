#include "tlab_error.h"

module FDM
    use TLab_Constants, only: wp, wi, roundoff_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_NONE
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLab_Grid, only: x, y, z, grid_dt
    use FDM_Derivative
    use FDM_Interpolate
    use FDM_Integral
    implicit none
    private

    type, public :: fdm_dt
        sequence
        character*8 name
        integer(wi) size
        logical :: uniform = .false.
        logical :: periodic = .false.
        real(wp) scale
        !
        real(wp), allocatable :: nodes(:)
        real(wp), allocatable :: jac(:, :)      ! grid spacing, Jacobian of 1. order derivative; need aux space for 2. order derivative
        !
        type(fdm_derivative_dt) :: der1
        type(fdm_derivative_dt) :: der2
        type(fdm_interpol_dt) :: intl

    end type fdm_dt

    type(fdm_dt), public, protected :: g(3)                    ! fdm derivative plans along 3 directions
    type(fdm_integral_dt), public, protected :: fdm_Int0(2)    ! fdm integral plans along Oy (ode for lambda = 0)

    public :: FDM_Initialize
    public :: FDM_CreatePlan

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Initialize(inifile)
        character(len=*), optional, intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer ig

        !########################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Main'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#SpaceOrder1=<CompactJacobian4/CompactJacobian6/CompactJacobian6Penta/CompactDirect4/CompactDirect6>')
        call TLab_Write_ASCII(bakfile, '#SpaceOrder2=<CompactJacobian4/CompactJacobian6/CompactJacobian6Hyper/CompactDirect4/CompactDirect6>')

        call ScanFile_Char(bakfile, inifile, block, 'SpaceOrder1', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, block, 'SpaceOrder', 'compactjacobian6', sRes) ! backwards compatibility
        if (trim(adjustl(sRes)) == 'compactjacobian4') then; g(1:3)%der1%mode_fdm = FDM_COM4_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6') then; g(1:3)%der1%mode_fdm = FDM_COM6_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6penta') then; g(1:3)%der1%mode_fdm = FDM_COM6_JACOBIAN_PENTA; 
        elseif (trim(adjustl(sRes)) == 'compactdirect4') then; g(1:3)%der1%mode_fdm = FDM_COM4_DIRECT; 
        elseif (trim(adjustl(sRes)) == 'compactdirect6') then; g(1:3)%der1%mode_fdm = FDM_COM6_DIRECT; 
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong SpaceOrder1 option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'SpaceOrder2', 'CompactJacobian6Hyper', sRes)
        if (trim(adjustl(sRes)) == 'compactjacobian4') then; g(1:3)%der2%mode_fdm = FDM_COM4_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6') then; g(1:3)%der2%mode_fdm = FDM_COM6_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6hyper') then; g(1:3)%der2%mode_fdm = FDM_COM6_JACOBIAN_HYPER; 
        elseif (trim(adjustl(sRes)) == 'compactdirect4') then; g(1:3)%der2%mode_fdm = FDM_COM4_DIRECT; 
        elseif (trim(adjustl(sRes)) == 'compactdirect6') then; g(1:3)%der2%mode_fdm = FDM_COM6_DIRECT; 
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong SpaceOrder2 option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (g(1)%der1%mode_fdm == FDM_COM6_JACOBIAN_PENTA) then     ! CFL_max depends on max[g(ig)%der1%mwn(:)]
            call TLab_Write_ASCII(wfile, __FILE__//'. CompactJacobian6Penta requires adjusted CFL-number depending on alpha and beta values.')
        end if

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[Grid]')
        call TLab_Write_ASCII(bakfile, '#XUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#YUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#ZUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#XPeriodic=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#YPeriodic=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#ZPeriodic=<yes/no>')

        g(1)%name = 'x'
        g(2)%name = 'y'
        g(3)%name = 'z'

        do ig = 1, 3
            call ScanFile_Char(bakfile, inifile, 'Grid', g(ig)%name(1:1)//'Uniform', 'void', sRes)
            if (trim(adjustl(sRes)) == 'yes') then; g(ig)%uniform = .true.
            else if (trim(adjustl(sRes)) == 'no') then; g(ig)%uniform = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Uniform '//g(ig)%name(1:1)//' grid')
                call TLab_Stop(DNS_ERROR_UNIFORMX)
            end if

            call ScanFile_Char(bakfile, inifile, 'Grid', g(ig)%name(1:1)//'Periodic', 'void', sRes)
            if (trim(adjustl(sRes)) == 'yes') then; g(ig)%periodic = .true.
            else if (trim(adjustl(sRes)) == 'no') then; g(ig)%periodic = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Periodic '//g(ig)%name(1:1)//' grid')
                call TLab_Stop(DNS_ERROR_IBC)
            end if

            ! consistency check
            if (g(ig)%periodic .and. (.not. g(ig)%uniform)) then
                call TLab_Write_ASCII(efile, __FILE__//'. Grid must be uniform in periodic direction.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end do
        x%periodic = g(1)%periodic
        y%periodic = g(2)%periodic
        z%periodic = g(3)%periodic
        
        !########################################################################
        ! Initializing fdm plan for derivatives
        call FDM_CreatePlan(x, g(1))
        call FDM_CreatePlan(y, g(2))
        call FDM_CreatePlan(z, g(3))

        ! ###################################################################
        ! Initializing fdm plans for first-order integrals (cases lambda = 0.0_wp)
        call FDM_Int1_Initialize(y%nodes(:), g(2)%der1, 0.0_wp, BCS_MIN, fdm_Int0(BCS_MIN))
        call FDM_Int1_Initialize(y%nodes(:), g(2)%der1, 0.0_wp, BCS_MAX, fdm_Int0(BCS_MAX))

        return
    end subroutine FDM_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_CreatePlan(x, g, locScale)
        type(grid_dt) :: x
        ! real(wp), intent(in) :: nodes(:)                            ! positions of the grid nodes
        type(fdm_dt), intent(inout) :: g                            ! fdm plan for derivatives
        real(wp), intent(in), optional :: locScale                  ! for consistency check

        ! -------------------------------------------------------------------
        integer(wi) i, nx

        ! ###################################################################
        ! Consistency check
        ! ###################################################################
        if (g%periodic .and. g%der1%mode_fdm == FDM_COM4_DIRECT) g%der1%mode_fdm = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%der1%mode_fdm == FDM_COM6_DIRECT) g%der1%mode_fdm = FDM_COM6_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%der2%mode_fdm == FDM_COM4_DIRECT) g%der2%mode_fdm = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%der2%mode_fdm == FDM_COM6_DIRECT) g%der2%mode_fdm = FDM_COM6_JACOBIAN_HYPER  ! they are the same for uniform grids.

        g%size = x%size !size(nodes)

        if (g%size > 1) then
            g%scale = x%nodes(g%size) - x%nodes(1)
            if (g%periodic) g%scale = g%scale*(1.0_wp + 1.0_wp/real(g%size - 1, wp))
        else
            g%scale = 1.0_wp  ! to avoid conditionals and NaN in some of the calculations below
        end if

        if (present(locScale)) then
! print *, abs((scale_loc - g%scale)/scale_loc)
            if (abs((locScale - g%scale)/g%scale) > roundoff_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Unmatched domain scale.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        ! ###################################################################
        nx = g%size                     ! # of nodes, for clarity below

        if (allocated(g%nodes)) deallocate (g%nodes)
        if (allocated(g%jac)) deallocate (g%jac)
        allocate (g%nodes(nx))
        allocate (g%jac(nx, 1 + 2))     ! I need 2 aux array to calculate the Jacobian for 2. order derivative; to be fixed

        if (nx == 1) then
            g%jac(:, :) = 1.0_wp
            return
        end if

        ! ###################################################################
        ! first-order derivative
        ! ###################################################################
        ! uniform grid to calculate Jacobian (used for the stencils below and also as grid spacing in the code).
        g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
        g%jac(:, 1) = 1.0_wp

        call FDM_Der1_Initialize(g%nodes, g%jac(:, 1), g%der1, periodic=.false., bcs_cases=[BCS_DD])

        ! Calculating derivative dxds into g%jac(:, 1)
        g%der1%periodic = .false.
        call FDM_Der1_Solve(1, BCS_NONE, g%der1, g%der1%lu, x%nodes, g%jac(:, 1), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        g%nodes(:) = x%nodes(1:nx)

        call FDM_Der1_Initialize(g%nodes, g%jac(:, 1), g%der1, periodic=g%periodic, bcs_cases=[BCS_DD, BCS_ND, BCS_DN, BCS_NN])

        if (g%periodic) g%der1%mwn(:) = g%der1%mwn(:)/g%jac(1, 1)           ! normalized by dx

        ! ###################################################################
        ! second-order derivative
        ! ###################################################################
        ! -------------------------------------------------------------------
        ! uniform grid to calculate Jacobian (needed to set up the Jacobian formulations).
        g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
        g%jac(:, 2) = 1.0_wp
        g%jac(:, 3) = 0.0_wp

        call FDM_Der2_Initialize(g%nodes, g%jac(:, 2:), g%der2, periodic=.false., uniform=.true.)

        ! Calculating derivative d2xds2 into g%jac(:, 3)
        g%der2%periodic = .false.
        call FDM_Der2_Solve(1, g%der2, g%der2%lu, x%nodes, g%jac(:, 3), g%jac(:, 2), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        g%nodes(:) = x%nodes(1:nx)
        g%jac(:, 2) = g%jac(:, 1)

        call FDM_Der2_Initialize(g%nodes, g%jac(:, 2:), g%der2, g%periodic, g%uniform)

        if (g%der2%periodic) g%der2%mwn(:) = g%der2%mwn(:)/(g%jac(1, 1)**2)      ! normalized by dx

        ! ###################################################################
        ! interpolation for staggered cases
        ! ###################################################################
        if (stagger_on) then
            if (g%periodic) then

                call FDM_Interpol_Initialize(x%nodes(:), g%jac(:, 1), g%intl, g%der1%mwn(:))

                ! else
                !     call TLab_Write_ASCII(efile, 'Staggered grid only along periodic directions.')
                !     call TLab_Stop(DNS_ERROR_UNDEVELOP)

            end if

        end if

        return
    end subroutine FDM_CreatePlan

end module FDM
