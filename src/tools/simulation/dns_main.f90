#include "tlab_error.h"

program DNS
    use TLab_Constants, only: ifile, efile, wfile, lfile, gfile, tag_flow, tag_scal, tag_part, tag_traj
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_WorkFlow, only: scal_on, flow_on
    use TLab_Time, only: itime, rtime
    use TLab_Arrays
    use TLab_Memory, only: imax, jmax, kmax, inb_scal, inb_flow, isize_field
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    ! use Thermodynamics, only: Thermodynamics_Initialize_Parameters
    use NavierStokes, only: NavierStokes_Initialize_Parameters, DNS_EQNS_ANELASTIC, DNS_EQNS_INCOMPRESSIBLE
    use NavierStokes, only: visc
    ! use Gravity, only: Gravity_Initialize
    ! use Rotation, only: Rotation_Initialize
    ! use Rotation, only: Rotation_Initialize
    ! use Radiation, only: Radiation_Initialize
    ! use Microphysics, only: Microphysics_Initialize
    ! use Chemistry, only: Chemistry_Initialize
    ! use SpecialForcing, only: SpecialForcing_Initialize
    ! use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use Tlab_Background, only: TLab_Initialize_Background!, pbg, rbg
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use OPR_Elliptic, only: OPR_Elliptic_Initialize
    use OPR_Burgers, only: OPR_Burgers_Initialize
    ! use OPR_FILTERS
    ! use PARTICLE_VARS
    ! use PARTICLE_ARRAYS
    ! use PARTICLE_PROCS
    use DNS_ARRAYS
    use DNS_LOCAL
    use DNS_Control
    use TimeMarching!, only: dtime
    ! use DNS_TOWER
    ! use IBM_VARS
    ! use PLANES
    ! use BOUNDARY_INFLOW
    ! use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    ! use DNS_STATISTICS, only: DNS_STATISTICS_INITIALIZE, DNS_STATISTICS_SPATIAL, DNS_STATISTICS_TEMPORAL, mean_flow, mean_scal
    ! use ParticleTrajectories
    ! use AVG_SCAL_ZT
    ! use AVG_PHASE
    ! use Avg_Spatial, only: IO_READ_AVG_SPATIAL, IO_WRITE_AVG_SPATIAL
    implicit none

    ! -------------------------------------------------------------------
    character(len=32) fname, str
    ! integer ig
    ! integer, parameter :: i0 = 0, i1 = 1
    real(wp) params(2)

    ! ###################################################################
    call system_clock(start_clock)
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif
    ! call Particle_Initialize_Parameters(ifile)
    ! call IBM_READ_INI(ifile)
    ! if (imode_ibm == 1) then
    !     call IBM_READ_CONSISTENCY_CHECK()
    ! end if

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    ! call Thermodynamics_Initialize_Parameters(ifile)
    ! call Gravity_Initialize(ifile)
    ! call Rotation_Initialize(ifile)
    ! call Radiation_Initialize(ifile)
    ! call Microphysics_Initialize(ifile)
    ! call LargeScaleForcing_Initialize(ifile)
    ! call Chemistry_Initialize(ifile)

    call TLab_Consistency_Check()

    call DNS_Initialize_Parameters(ifile)

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    ! call SpecialForcing_Initialize(ifile)

    call TLab_Initialize_Background(ifile)

    call TLab_Allocate_Real(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
    call TLab_Allocate_Real(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

    ! call ParticleTrajectories_Initialize(ifile)
    ! call Particle_Initialize_Memory(__FILE__)
    ! call TLab_Allocate_Real(__FILE__, l_hq, [isize_part, inb_part], 'part-rhs')

    ! call DNS_STATISTICS_INITIALIZE()

    ! call PLANES_INITIALIZE()

    ! if (PhAvg%active) then
    !     call AvgPhaseInitializeMemory(__FILE__, nitera_save)
    ! end if

    ! if (use_tower) then
    !     call DNS_TOWER_INITIALIZE(tower_stride)
    ! end if

    ! if (imode_ibm == 1) then
    !     call IBM_ALLOCATE(__FILE__)
    ! end if

    ! ###################################################################
    ! Initialize operators
    ! ###################################################################
    call OPR_Burgers_Initialize(ifile)

    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    ! call OPR_Filter_Initialize_Parameters(ifile)
    ! do ig = 1, 3
    !     call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
    !     call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    ! end do

    ! ###################################################################
    ! Initialize fields
    ! ###################################################################
    itime = nitera_first

    visc_stop = visc ! Value read in ifile

    if (scal_on) then
        write (fname, *) nitera_first; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
    end if

    write (fname, *) nitera_first; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
    call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:2))
    rtime = params(1); visc = params(2)

    ! call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)  ! Initialize diagnostic thermodynamic quantities

    ! if (part%type /= PART_TYPE_NONE) then
    !     write (fname, *) nitera_first; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
    !     call IO_READ_PARTICLE(fname, l_g, l_q)
    !     call Particle_Initialize_Fields()
    ! end if

    ! ###################################################################
    ! Initialize change in viscosity
    ! ###################################################################
    flag_viscosity = .false.
    if (visc /= visc_stop) then
        write (str, *) visc
        call TLab_Write_ASCII(lfile, 'Changing original viscosity '//trim(adjustl(str))//' to new value.')
        if (visc_time > 0.0_wp) then
            visc_rate = (visc_stop - visc)/visc_time
            visc_time = rtime + visc_time                 ! Stop when this time is reached
            flag_viscosity = .true.
        else
            visc = visc_stop
        end if
    end if

    ! ###################################################################
    ! Initialize data for boundary conditions
    ! ###################################################################
    ! call BOUNDARY_BUFFER_INITIALIZE(q, s, txc)

    call BOUNDARY_BCS_INITIALIZE()

    ! ! ###################################################################
    ! ! Initialize IBM
    ! ! ###################################################################
    ! if (imode_ibm == 1) then
    !     call IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
    !     call IBM_BCS_FIELD_COMBINED(i0, q)
    !     if (scal_on) call IBM_INITIALIZE_SCAL(i1, s)
    ! end if

    ! ###################################################################
    ! Check
    ! ###################################################################

    ! ###################################################################
    ! Initialize time marching scheme
    ! ###################################################################
    call TMarch_Initialize(ifile)
    call TMarch_Courant()

    ! ###################################################################
    ! Check-pointing: Initialize logfiles, write header & first line
    ! ###################################################################
    call DNS_Control_Initialize(ifile)

    ! ###################################################################
    ! Do simulation: Integrate equations
    ! ###################################################################
    itime = nitera_first

    write (str, *) itime
    call TLab_Write_ASCII(lfile, 'Starting time integration at It'//trim(adjustl(str))//'.')

    do
        if (itime >= nitera_last) exit
        if (int(logs_data(1)) /= 0) exit
        call TMarch_RungeKutta()
        itime = itime + 1
        rtime = rtime + dtime
        ! if (mod(itime - nitera_first, nitera_filter) == 0) then
        !     call DNS_FILTER()
        !     if (imode_ibm == 1) then
        !         call IBM_BCS_FIELD_COMBINED(i0, q) ! apply IBM BCs
        !         if (scal_on) call IBM_INITIALIZE_SCAL(i0, s)
        !     end if
        ! end if

        if (flag_viscosity) then                ! Change viscosity if necessary
            visc = visc + visc_rate*dtime
            if (rtime > visc_time) then
                visc = visc_stop                ! Fix new value without any roundoff
                flag_viscosity = .false.
            end if
        end if

        call TMarch_Courant()

        ! -------------------------------------------------------------------
        ! The rest: Logging, postprocessing and check-pointing
        ! -------------------------------------------------------------------
        call DNS_BOUNDS_CONTROL()
        call DNS_OBS_CONTROL()
        if (mod(itime - nitera_first, nitera_log) == 0 .or. int(logs_data(1)) /= 0) then
            call DNS_LOGS()
            if (dns_obs_log /= OBS_TYPE_NONE) then
                call DNS_OBS()
            end if
        end if

        ! if (PhAvg%active) then
        !     if (mod(itime, PhAvg%stride) == 0) then
        !         call AvgPhaseSpace(wrk2d, inb_flow, itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, 1)
        !         call AvgPhaseSpace(wrk2d, inb_scal, itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, 2)
        !         ! Pressure is taken from the RHS subroutine
        !         ! call AvgPhaseSpace(wrk2d, 6       , itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, 8)
        !         call AvgPhaseStress(q, itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride)
        !         if (mod(itime - nitera_first, nitera_save) == 0) then
        !             call IO_Write_AvgPhase(avg_planes, inb_flow, IO_FLOW, nitera_save, PhAvg%stride, avgu_name, 1, avg_flow)
        !             call IO_Write_AvgPhase(avg_planes, inb_scal, IO_SCAL, nitera_save, PhAvg%stride, avgs_name, 2, avg_scal)
        !             call IO_Write_AvgPhase(avg_planes, 1, IO_SCAL, nitera_save, PhAvg%stride, avgp_name, 4, avg_p)
        !             call IO_Write_AvgPhase(avg_planes, 6, IO_FLOW, nitera_save, PhAvg%stride, avgstr_name, 8, avg_stress)

        !             call AvgPhaseResetVariable()
        !         end if
        !     end if
        ! end if

        ! if (use_tower) then
        !     call DNS_TOWER_ACCUMULATE(q, 1, wrk1d)
        !     call DNS_TOWER_ACCUMULATE(s, 2, wrk1d)
        ! end if
        ! if (imode_traj /= TRAJ_TYPE_NONE) then
        !     call ParticleTrajectories_Accumulate()
        ! end if

        if (mod(itime - nitera_first, nitera_save) == 0 .or. &      ! Check-pointing: Save restart files
            itime == nitera_last .or. int(logs_data(1)) /= 0 .or. & ! Secure that one restart file is saved
            wall_time > nruntime_sec) then                          ! If max runtime of the code is reached

            if (flow_on) then
                write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
                io_header_q(1)%params(1) = rtime
                call IO_Write_Fields(fname, imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))
            end if

            if (scal_on) then
                write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
                io_header_s(:)%params(1) = rtime
                call IO_Write_Fields(fname, imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))
            end if

            ! if (use_tower) then
            !     call DNS_TOWER_WRITE(wrk3d)
            ! end if

            ! if (part%type /= PART_TYPE_NONE) then
            !     write (fname, *) itime; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
            !     call IO_WRITE_PARTICLE(fname, l_g, l_q)
            !     if (imode_traj /= TRAJ_TYPE_NONE) then
            !         write (fname, *) itime; fname = trim(adjustl(tag_traj))//trim(adjustl(fname))
            !         call ParticleTrajectories_Write(fname)
            !     end if
            ! end if

        end if

        ! if (mod(itime - nitera_first, nitera_pln) == 0) then
        !     call PLANES_SAVE()
        ! end if

        if (wall_time > nruntime_sec) then
            write (str, *) wall_time
            ! write to efile so that job is not resubmitted
            call TLab_Write_ASCII(efile, 'Maximum walltime of '//trim(adjustl(str))//' seconds is reached.')
            exit
        end if
    end do

    ! ###################################################################
    call TLab_Stop(int(logs_data(1)))

end program DNS
