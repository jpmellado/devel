#include "tlab_error.h"

module DNS_ARRAYS
    use TLab_Constants, only: wp
    implicit none

    real(wp), allocatable, target :: hq(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, target :: hs(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, target :: l_hq(:, :)     ! Right-hand sides Lagrangian fields

end module DNS_ARRAYS

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

    logical :: remove_divergence    ! Remove residual divergence every time step

    ! Variable viscosity
    logical :: flag_viscosity
    real(wp) :: visc_stop, visc_time, visc_rate

contains
    ! ###################################################################
    ! ###################################################################
    subroutine DNS_READ_LOCAL(inifile)
        use TLab_Constants, only: wp, wi, big_wp, efile, lfile, wfile
        use FDM, only: g
        use TLab_Memory, only: inb_scal, inb_txc
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        use NavierStokes, only: DNS_EQNS_ANELASTIC, DNS_EQNS_INCOMPRESSIBLE
        ! use BOUNDARY_BUFFER
        use BOUNDARY_BCS
        ! use PARTICLE_VARS
        ! use DNS_STATISTICS, only: stats_averages, stats_pdfs, stats_intermittency
        ! use PLANES
        ! use IBM_VARS, only: imode_ibm, ibm_geo
        ! use AVG_PHASE
        ! use OPR_Filters, only: FilterDomain, PressureFilter, DNS_FILTER_NONE
        ! use Discrete, only: Discrete_ReadBlock

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, eStr
        character(len=512) sRes
        integer is

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        ! ###################################################################
        block = 'Main'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#TermDivergence=<none/remove>')

        call ScanFile_Char(bakfile, inifile, 'Main', 'TermDivergence', 'remove', sRes)
        if (trim(adjustl(sRes)) == 'none') then; remove_divergence = .false.
        else if (trim(adjustl(sRes)) == 'remove') then; remove_divergence = .true.
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TermDivergence option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        block = 'Iteration'
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
        call TLab_Write_ASCII(bakfile, '#IteraLog=<value>')

        call ScanFile_Int(bakfile, inifile, 'Iteration', 'Start', '0', nitera_first)
        call ScanFile_Int(bakfile, inifile, 'Iteration', 'End', '0', nitera_last)
        call ScanFile_Int(bakfile, inifile, 'Iteration', 'Restart', '50', nitera_save)
        call ScanFile_Int(bakfile, inifile, 'Iteration', 'Statistics', '50', nitera_stats)
        call ScanFile_Int(bakfile, inifile, 'Iteration', 'Saveplanes', '-1', nitera_pln)
        call ScanFile_Real(bakfile, inifile, 'Iteration', 'Runtime', '10000000', nruntime_sec)
        call ScanFile_Int(bakfile, inifile, 'Iteration', 'IteraLog', '10', nitera_log)

        ! ! Domain Filter (Should we move it to Iteration?)
        ! call ScanFile_Int(bakfile, inifile, 'Filter', 'Step', '0', nitera_filter)
        ! if (nitera_filter == 0) FilterDomain(:)%type = DNS_FILTER_NONE

        ! ###################################################################
        block = 'BoundaryConditions'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#ScalarImin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarJmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarKmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarSfcTypeJmin=<static/linear>')
        call TLab_Write_ASCII(bakfile, '#ScalarSfcTypeJmax=<static/linear>')
        call TLab_Write_ASCII(bakfile, '#ScalarCouplingJmin=<value>')
        call TLab_Write_ASCII(bakfile, '#ScalarCouplingJmax=<value>')
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
                call TLab_Write_ASCII(efile, &
                                      'DNS_READ_LOCAL. Interactive BC at jmin not implemented for non-Dirichlet BC')
                call TLab_Stop(DNS_ERROR_JBC)
            end if
            if (BcsScalJmax%type(is) /= DNS_BCS_DIRICHLET .and. &
                BcsScalJmax%SfcType(is) /= DNS_SFC_STATIC) then
                write (*, *) BcsScalJmax%type(is), BcsScalJmax%SfcType(is), BcsScalJmax%cpl(is)
                call TLab_Write_ASCII(efile, &
                                      'DNS_READ_LOCAL. Interactive BC at jmax not implemented for non-Dirichlet BC')
                call TLab_Stop(DNS_ERROR_JBC)
            end if
        end do

        ! ###################################################################
        block = 'BufferZone'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/relaxation/filter/both>')
        call TLab_Write_ASCII(bakfile, '#LoadBuffer=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#PointsJmin=<value>')
        call TLab_Write_ASCII(bakfile, '#PointsJmax=<value>')
        call TLab_Write_ASCII(bakfile, '#PointsImin=<value>')
        call TLab_Write_ASCII(bakfile, '#PointsImax=<value>')
        call TLab_Write_ASCII(bakfile, '#ParametersJmin=<values>')
        call TLab_Write_ASCII(bakfile, '#ParametersJmax=<values>')
        call TLab_Write_ASCII(bakfile, '#ParametersImin=<values>')
        call TLab_Write_ASCII(bakfile, '#ParametersImax=<values>')
        call TLab_Write_ASCII(bakfile, '#HardValuesJmin=<values>')
        call TLab_Write_ASCII(bakfile, '#HardValuesJmax=<values>')
        call TLab_Write_ASCII(bakfile, '#HardValuesImin=<values>')
        call TLab_Write_ASCII(bakfile, '#HardValuesImax=<values>')

        !         call ScanFile_Char(bakfile, inifile, 'BufferZone', 'Type', 'none', sRes)
        !         if (trim(adjustl(sRes)) == 'none') then; BuffType = DNS_BUFFER_NONE
        !         else if (trim(adjustl(sRes)) == 'relaxation') then; BuffType = DNS_BUFFER_RELAX
        !         else if (trim(adjustl(sRes)) == 'filter') then; BuffType = DNS_BUFFER_FILTER
        !         else if (trim(adjustl(sRes)) == 'both') then; BuffType = DNS_BUFFER_BOTH
        !         else
        !             call TLab_Write_ASCII(efile, 'DNS_READ_LOCAL. Wrong BufferType option.')
        !             call TLab_Stop(DNS_ERROR_OPTION)
        !         end if

        ! ! Load buffer if used also by BCs
        !         call ScanFile_Char(bakfile, inifile, 'BufferZone', 'LoadBuffer', 'no', sRes)
        !         if (trim(adjustl(sRes)) == 'yes') then; BuffLoad = .true.
        !         else; BuffLoad = .false.
        !         end if

        ! ! Sizes; read always because allocation checks if # points is zero
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsUImin', '0', BuffFlowImin%size)
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsUImax', '0', BuffFlowImax%size)
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsUJmin', '0', BuffFlowJmin%size)
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsUJmax', '0', BuffFlowJmax%size)

        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsSImin', '0', BuffScalImin%size)
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsSImax', '0', BuffScalImax%size)
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsSJmin', '0', BuffScalJmin%size)
        !         call ScanFile_Int(bakfile, inifile, 'BufferZone', 'PointsSJmax', '0', BuffScalJmax%size)

        !         BuffScalImin%type = BuffType; BuffFlowImin%type = BuffType ! So far, all the same
        !         BuffScalImax%type = BuffType; BuffFlowImax%type = BuffType
        !         BuffScalJmin%type = BuffType; BuffFlowJmin%type = BuffType
        !         BuffScalJmax%type = BuffType; BuffFlowJmax%type = BuffType

        !         if (BuffScalImin%size /= BuffFlowImin%size .or. &
        !             BuffScalImax%size /= BuffFlowImax%size .or. &
        !             BuffScalJmin%size /= BuffFlowJmin%size .or. &
        !             BuffScalJmax%size /= BuffFlowJmax%size) then ! Because of io_subarray
        !             call TLab_Write_ASCII(wfile, 'DNS_READ_LOCAL. Buffer zone sizes must be equal in flow and scal.')
        !             call TLab_Stop(DNS_ERROR_OPTION)
        !         end if

        ! ! Parameters
        !         if (BuffType /= DNS_BUFFER_NONE) then
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UImin', BuffFlowImin, inb_flow)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UImax', BuffFlowImax, inb_flow)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UJmin', BuffFlowJmin, inb_flow)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'UJmax', BuffFlowJmax, inb_flow)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SImin', BuffScalImin, inb_scal)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SImax', BuffScalImax, inb_scal)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SJmin', BuffScalJmin, inb_scal)
        !             call BOUNDARY_BUFFER_READBLOCK(bakfile, inifile, 'SJmax', BuffScalJmax, inb_scal)
        !         end if

        ! ! Make sure there is array space for reference mean drift
        !         if (BcsDrift) then
        !             BuffFlowJmin%size = max(BuffFlowJmin%size, 1); BuffFlowJmax%size = max(BuffFlowJmax%size, 1)
        !             BuffScalJmin%size = BuffFlowJmin%size; BuffScalJmax%size = BuffFlowJmax%size
        !             BuffScalImin%size = BuffFlowImin%size; BuffScalImax%size = BuffFlowImax%size
        !         end if

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
        !             call TLab_Write_ASCII(wfile, 'DNS_READ_LOCAL. Cannot read stride for towers; set to 0,0,0.')
        !         end if

        ! ###################################################################
        ! Statistics Control
        ! ###################################################################
        ! call TLab_Write_ASCII(bakfile, '#')
        ! call TLab_Write_ASCII(bakfile, '#[Statsitics]')
        ! call TLab_Write_ASCII(bakfile, '#Averages=<yes/no>')
        ! call TLab_Write_ASCII(bakfile, '#Pdfs=<yes/no>')
        ! call TLab_Write_ASCII(bakfile, '#Intermittency=<yes/no>')

        ! call ScanFile_Char(bakfile, inifile, 'Statistics', 'Averages', 'yes', sRes)
        ! if (trim(adjustl(sRes)) == 'yes') then; stats_averages = .true.
        ! else; stats_averages = .false.; end if

        ! call ScanFile_Char(bakfile, inifile, 'Statistics', 'Pdfs', 'yes', sRes)
        ! if (trim(adjustl(sRes)) == 'yes') then; stats_pdfs = .true.
        ! else; stats_pdfs = .false.; end if

        ! call ScanFile_Char(bakfile, inifile, 'Statistics', 'Intermittency', 'yes', sRes)
        ! if (trim(adjustl(sRes)) == 'yes') then; stats_intermittency = .true.
        ! else; stats_intermittency = .false.; end if

        ! ! ###################################################################
        ! ! Phase Averaging
        ! ! ###################################################################
        ! call ScanFile_Int(bakfile, inifile, 'Iteration', 'PhaseAvg', '0', PhAvg%stride)
        ! if (PhAvg%stride > 0) PhAvg%active = .true.

        ! ###################################################################
        ! Final initialization and consistency check
        ! ###################################################################
        if (nitera_first > nitera_last) then
            call TLab_Write_ASCII(efile, 'DNS_READ_LOCAL. Not started because nitera_first > nitera_last.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

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
        !                 call TLab_Write_ASCII(efile, 'DNS_READ_LOCAL. Not enough space in wrk3d array to handle tower information. Increase strides.')
        !                 call TLab_Stop(DNS_ERROR_UNDEVELOP)
        !             end if
        !         end if

        ! ! -------------------------------------------------------------------
        ! ! Pressure staggering and filtering
        ! ! -------------------------------------------------------------------
        !         if (stagger_on .or. any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        !             if (.not. (imode_rhs == EQNS_RHS_COMBINED)) then
        !                 call TLab_Write_ASCII(efile, 'DNS_READ_LOCAL. Horizontal pressure staggering or Pressure filter not implemented for this RHS type.')
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
    end subroutine DNS_READ_LOCAL

end module DNS_LOCAL
