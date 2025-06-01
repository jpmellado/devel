#include "tlab_error.h"

module DNS_Arrays
    use TLab_Constants, only: wp
    implicit none

    real(wp), allocatable, target :: hq(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, target :: hs(:, :)       ! Right-hand sides Eulerian fields
    ! real(wp), allocatable, target :: l_hq(:, :)     ! Right-hand sides Lagrangian fields

    real(wp), pointer :: p_hq(:, :, :, :) => null()
    real(wp), pointer :: p_hs(:, :, :, :) => null()

end module DNS_Arrays

! ###################################################################
! ###################################################################
module DNS_LOCAL
    use TLab_Constants, only: wp, wi, sp
    ! use TLab_Constants, only: MAX_PATH_LENGTH
    implicit none

    integer :: nitera_first     ! First iteration in current run
    integer :: nitera_last      ! Last iteration in current run
    integer :: nitera_save      ! Iteration step to check-point: save restart files
    integer :: nitera_stats     ! Iteration step to check-point: save statistical data
    integer :: nitera_pln       ! Iteration step to save planes
    integer :: nitera_filter    ! Iteration step for domain filter, if any
    integer :: nitera_log           ! Iteration step for data logger with simulation information

    real(wp) :: nruntime_sec     ! Maximum runtime of the simulation in seconds
    real(wp) :: wall_time        ! Actual elapsed time during the simulation in seconds
    integer :: start_clock      ! Starting time of the simulation on the system

    ! Variable viscosity
    logical :: flag_viscosity
    real(wp) :: visc_stop, visc_time, visc_rate

contains
    ! ###################################################################
    ! ###################################################################
    subroutine DNS_Initialize_Parameters(inifile)
        use TLab_Constants, only: wp, wi, big_wp, efile, lfile, wfile
        use FDM, only: g
        use TLab_Memory, only: inb_scal, inb_txc
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        use NavierStokes, only: DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
        ! use BOUNDARY_BUFFER
        use BOUNDARY_BCS
        ! use PARTICLE_VARS
        ! use PLANES
        ! use OPR_Filters, only: FilterDomain, PressureFilter, DNS_FILTER_NONE
        ! use Discrete, only: Discrete_ReadBlock

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        ! character(len=512) sRes
        integer is

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Start=<integral start time>')
        call TLab_Write_ASCII(bakfile, '#End=<integral stop time>')
        call TLab_Write_ASCII(bakfile, '#Restart=<restart time step>')
        call TLab_Write_ASCII(bakfile, '#Statistics=<statistics time step>')
        call TLab_Write_ASCII(bakfile, '#Saveplanes=<value>')
        call TLab_Write_ASCII(bakfile, '#RunAvera=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#Runtime=<seconds>')
        call TLab_Write_ASCII(bakfile, '#Logs=<value>')

        call ScanFile_Int(bakfile, inifile, block, 'Start', '0', nitera_first)
        call ScanFile_Int(bakfile, inifile, block, 'End', '0', nitera_last)
        call ScanFile_Int(bakfile, inifile, block, 'Restart', '50', nitera_save)
        call ScanFile_Int(bakfile, inifile, block, 'Statistics', '50', nitera_stats)
        call ScanFile_Int(bakfile, inifile, block, 'Saveplanes', '-1', nitera_pln)
        call ScanFile_Int(bakfile, inifile, block, 'Logs', '10', nitera_log)
        call ScanFile_Real(bakfile, inifile, block, 'Runtime', '10000000', nruntime_sec)

        ! ! Domain Filter (Should we move it to Iteration?)
        ! call ScanFile_Int(bakfile, inifile, 'Filter', 'Step', '0', nitera_filter)
        ! if (nitera_filter == 0) FilterDomain(:)%type = DNS_FILTER_NONE

        if (nitera_first > nitera_last) then
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Not started because nitera_first > nitera_last.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        block = 'BoundaryConditions'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#ScalarImin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarJmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarKmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarSfcTypeKmin=<static/linear>')
        call TLab_Write_ASCII(bakfile, '#ScalarSfcTypeKmax=<static/linear>')
        call TLab_Write_ASCII(bakfile, '#ScalarCouplingKmin=<value>')
        call TLab_Write_ASCII(bakfile, '#ScalarCouplingKmax=<value>')
        call TLab_Write_ASCII(bakfile, '#VelocityImin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#VelocityJmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#VelocityKmin=<none/dirichlet/neumman>')

        ! -------------------------------------------------------------------
        ! Scalar terms (including surface model at vertical boundaries)
        ! -------------------------------------------------------------------
        BcsScalImin%type(:) = DNS_BCS_NONE; BcsScalImax%type(:) = DNS_BCS_NONE
        if (.not. g(1)%periodic) then
            call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Imin', BcsScalImin)
            call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Imax', BcsScalImax)
        end if

        BcsScalJmin%type(:) = DNS_BCS_NONE; BcsScalJmax%type(:) = DNS_BCS_NONE
        if (.not. g(2)%periodic) then
            call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Jmin', BcsScalJmin)
            call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Jmax', BcsScalJmax)
        end if

        BcsScalKmin%type(:) = DNS_BCS_NONE; BcsScalKmax%type(:) = DNS_BCS_NONE
        if (.not. g(3)%periodic) then
            call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Kmin', BcsScalKmin)
            call BOUNDARY_BCS_SCAL_READBLOCK(bakfile, inifile, 'Kmax', BcsScalKmax)
        end if

        ! -------------------------------------------------------------------
        ! Velocity terms / Euler part in compressible mode
        ! -------------------------------------------------------------------
        call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Imin', BcsFlowImin)
        call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Imax', BcsFlowImax)
        call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Jmin', BcsFlowJmin)
        call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Jmax', BcsFlowJmax)
        call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Kmin', BcsFlowKmin)
        call BOUNDARY_BCS_FLOW_READBLOCK(bakfile, inifile, 'Kmax', BcsFlowKmax)

        ! -------------------------------------------------------------------
        ! Boundary conditions
        ! -------------------------------------------------------------------
        ! Make sure periodic BCs are not modified
        if (g(1)%periodic) then; 
            BcsFlowImin%type(:) = DNS_BCS_NONE; BcsFlowImax%type(:) = DNS_BCS_NONE
            BcsScalImin%type(:) = DNS_BCS_NONE; BcsScalImax%type(:) = DNS_BCS_NONE
        end if
        if (g(2)%periodic) then; 
            BcsFlowJmin%type(:) = DNS_BCS_NONE; BcsFlowJmax%type(:) = DNS_BCS_NONE
            BcsScalJmin%type(:) = DNS_BCS_NONE; BcsScalJmax%type(:) = DNS_BCS_NONE
        end if
        if (g(3)%periodic) then; 
            BcsFlowKmin%type(:) = DNS_BCS_NONE; BcsFlowKmax%type(:) = DNS_BCS_NONE
            BcsScalKmin%type(:) = DNS_BCS_NONE; BcsScalKmax%type(:) = DNS_BCS_NONE
        end if

        ! -------------------------------------------------------------------
        ! Interactive Boundary conditions
        ! -------------------------------------------------------------------
        do is = 1, inb_scal
            if (BcsScalJmin%type(is) /= DNS_BCS_DIRICHLET .and. &
                BcsScalJmin%SfcType(is) /= DNS_SFC_STATIC) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Interactive BC at kmin not implemented for non-Dirichlet BC')
                call TLab_Stop(DNS_ERROR_JBC)
            end if
            if (BcsScalJmax%type(is) /= DNS_BCS_DIRICHLET .and. &
                BcsScalJmax%SfcType(is) /= DNS_SFC_STATIC) then
                write (*, *) BcsScalJmax%type(is), BcsScalJmax%SfcType(is), BcsScalJmax%cpl(is)
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Interactive BC at kmax not implemented for non-Dirichlet BC')
                call TLab_Stop(DNS_ERROR_JBC)
            end if
        end do

        ! ###################################################################
        ! Viscosity Control
        ! ###################################################################
        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[ViscChange]')
        call TLab_Write_ASCII(bakfile, '#Time=<time>')

        call ScanFile_Real(bakfile, inifile, 'ViscChange', 'Time', '0.0', visc_time)

        ! ! ###################################################################
        ! ! Save planes to disk
        ! ! ###################################################################
        !         call TLab_Write_ASCII(bakfile, '#')
        !         call TLab_Write_ASCII(bakfile, '#[SavePlanes]')

        !         call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesI', iplanes)
        !         call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesJ', jplanes)
        !         call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesK', kplanes)

        ! ! ###################################################################
        ! ! Save lines to disk
        ! ! ###################################################################
        !         call TLab_Write_ASCII(bakfile, '#')
        !         call TLab_Write_ASCII(bakfile, '#[SaveTowers]')
        !         call TLab_Write_ASCII(bakfile, 'Stride=<value_i,value_j,value_k>')

        !         call ScanFile_Char(bakfile, inifile, 'SaveTowers', 'Stride', '0,0,0', sRes)
        !         idummy = 3; call LIST_INTEGER(sRes, idummy, tower_stride)
        !         if (idummy /= 3) then
        !             tower_stride(:) = 0
        !             call TLab_Write_ASCII(bakfile, 'Stride=0,0,0')
        !             call TLab_Write_ASCII(wfile, 'DNS_Initialize_Parameters. Cannot read stride for towers; set to 0,0,0.')
        !         end if

        ! ###################################################################
        ! Final initialization and consistency check
        ! ###################################################################
        ! Avoid dividing by zero in time_integration routine
        if (nitera_save <= 0) nitera_save = nitera_last - nitera_first + 1
        if (nitera_stats <= 0) nitera_stats = nitera_last - nitera_first + 1
        if (nitera_log <= 0) nitera_log = nitera_last - nitera_first + 1
        if (nitera_pln <= 0) nitera_pln = nitera_last - nitera_first + 1
        if (nitera_filter <= 0) nitera_filter = nitera_last - nitera_first + 1

        ! ! -------------------------------------------------------------------
        ! ! Towers information
        ! ! So far, the use or not use of tower information is
        ! ! set by the stride information.
        ! ! -------------------------------------------------------------------
        !         use_tower = .false.
        !         if (minval(tower_stride) > 0) then; use_tower = .true.; end if

        ! ! We need space
        !         if (use_tower) then
        !             idummy = tower_stride(1)*tower_stride(2)*tower_stride(3)
        !             if (idummy < 5) then
        !                 call TLab_Write_ASCII(efile, 'DNS_Initialize_Parameters. Not enough space in wrk3d array to handle tower information. Increase strides.')
        !                 call TLab_Stop(DNS_ERROR_UNDEVELOP)
        !             end if
        !         end if

        ! ! -------------------------------------------------------------------
        ! ! Pressure staggering and filtering
        ! ! -------------------------------------------------------------------
        !         if (stagger_on .or. any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        !             if (.not. (imode_rhs == EQNS_RHS_COMBINED)) then
        !                 call TLab_Write_ASCII(efile, 'DNS_Initialize_Parameters. Horizontal pressure staggering or Pressure filter not implemented for this RHS type.')
        !                 call TLab_Stop(DNS_ERROR_UNDEVELOP)
        !             end if
        !         end if

        ! ! -------------------------------------------------------------------
        ! ! IBM Data
        ! ! -------------------------------------------------------------------
        !         if (imode_ibm == 1) then
        !             if (.not. (imode_rhs == EQNS_RHS_COMBINED)) then
        !                 call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. IBM only implemented for combined rhs mode.')
        !                 call TLab_Stop(DNS_ERROR_UNDEVELOP)
        !             end if

        !             do is = 1, inb_scal
        !                 if ((BcsScalJmin%type(is) /= DNS_BCS_DIRICHLET) .or. (BcsScalJmin%SfcType(is) /= DNS_SFC_STATIC) .or. &
        !                     (ibm_geo%mirrored .and. ((BcsScalJmax%type(is) /= DNS_BCS_DIRICHLET) .or. (BcsScalJmax%SfcType(is) /= DNS_SFC_STATIC)))) then
        !                     call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Wrong scalar BCs.')
        !                     call TLab_Stop(DNS_ERROR_UNDEVELOP)
        !                 end if
        !             end do
        !             do is = 1, 3
        !                 if ((BcsFlowJmin%type(is) /= DNS_BCS_DIRICHLET) .or. &
        !                     (ibm_geo%mirrored .and. (BcsFlowJmax%type(is) /= DNS_BCS_DIRICHLET))) then
        !                     call TLab_Write_ASCII(efile, 'IBM_READ_INI. IBM. Wrong Flow BCs.')
        !                     call TLab_Stop(DNS_ERROR_UNDEVELOP)
        !                 end if
        !             end do

        !         end if

        ! -------------------------------------------------------------------
        ! Array sizes
        ! -------------------------------------------------------------------

        inb_txc = 9

        ! ! isize_wrk3d = max(isize_wrk3d, g_inf(1)%size*g_inf(2)%size*g_inf(3)%size)
        !         if (use_tower) then
        !             isize_wrk3d = max(isize_wrk3d, nitera_save*(g(2)%size + 2))
        !         end if

        return
    end subroutine DNS_Initialize_Parameters

end module DNS_LOCAL
