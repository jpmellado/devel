#include "tlab_error.h"

program DNS
    use TLab_Constants, only: ifile, efile, wfile, lfile, gfile, tag_flow, tag_scal
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_WorkFlow, only: scal_on, flow_on
    use TLab_Time, only: itime, rtime
    use TLab_Arrays
    use TLab_Memory, only: imax, jmax, kmax, inb_scal, inb_flow
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    use NavierStokes, only: NavierStokes_Initialize_Parameters, DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
    use Thermodynamics, only: Thermo_Initialize
    use NavierStokes, only: visc
    use Gravity, only: Gravity_Initialize
    use SpecialForcing, only: SpecialForcing_Initialize
    ! use Rotation, only: Rotation_Initialize
    use Microphysics, only: Microphysics_Initialize
    use Radiation, only: Radiation_Initialize
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use OPR_Partial, only: OPR_Partial_Initialize
    use Tlab_Background, only: TLab_Initialize_Background!, pbg, rbg
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use OPR_Elliptic, only: OPR_Elliptic_Initialize
    use NSE_Burgers, only: NSE_Burgers_Initialize
    ! use OPR_FILTERS
    ! use PARTICLE_VARS
    ! use PARTICLE_ARRAYS
    ! use PARTICLE_PROCS
    use DNS_LOCAL
    use DNS_Control
    use TimeMarching
    ! use DNS_TOWER
    ! use PLANES
    use BOUNDARY_BCS
    use Buffer, only: Buffer_Initialize
    use Statistics, only: Statistics_Initialize, Statistics_Compute
    ! use ParticleTrajectories
    implicit none

    ! -------------------------------------------------------------------
    character(len=32) fname, str
    real(wp) params(2)

    ! ###################################################################
    call system_clock(start_clock)
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermo_Initialize(ifile)

    call Gravity_Initialize(ifile)
    call SpecialForcing_Initialize(ifile)
    ! call Rotation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)

    call TLab_Consistency_Check()

    call DNS_Initialize_Parameters(ifile)

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call OPR_Partial_Initialize()
    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    call TLab_Initialize_Background(ifile)
    call NSE_Burgers_Initialize(ifile)

    call Statistics_Initialize(ifile)

    ! call PLANES_INITIALIZE()


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

    call TLab_Diagnostic(imax, jmax, kmax, s)  ! Initialize diagnostic thermodynamic quantities

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
    call Buffer_Initialize(ifile)!q, s, txc)

    call BOUNDARY_BCS_INITIALIZE()

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

        if (mod(itime - nitera_first, nitera_stats) == 0) then      ! Calculate statistics
            call Statistics_Compute()
        end if

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
