#include "tlab_error.h"

module TLab_WorkFlow
    use TLab_Constants, only: sp, wp, wi, longi, lfile, efile, fmt_r
#ifdef USE_OPENMP
    use OMP_LIB
#endif
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro, ims_npro
    use TLabMPI_VARS, only: ims_time_max, ims_time_min, ims_time_trans
    use TLabMPI_VARS, only: ims_err
#endif
    implicit none
    private
    save

    character*128 :: line

    integer, public :: imode_verbosity = 1      ! level of verbosity used in log files

    logical, public :: flow_on = .true.         ! calculate flow parts of the code
    logical, public :: scal_on = .true.         ! calculate scal parts of the code
    logical, public :: stagger_on = .false.     ! horizontal staggering of pressure

    public :: TLab_Start
    public :: TLab_Stop
    public :: TLab_Write_ASCII

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Start()
        use TLab_OpenMP

        character*10 clock(2)

        !#####################################################################
        ! Inititalize MPI parallel mode
#ifdef USE_MPI
        call MPI_INIT(ims_err)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
        call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)

        write (line, *) ims_npro
        line = 'Number of MPI tasks '//trim(adjustl(line))
        call TLab_Write_ASCII(lfile, line)

        if (ims_npro == 0) then
            call TLab_Write_ASCII(efile, __FILE__//'. Number of processors is zero.')
            call TLab_Stop(DNS_ERROR_MINPROC)
        end if

        ims_time_min = MPI_WTIME()
        ims_time_trans = 0.0_wp

#endif
        !########################################################################
        ! First output
        call date_and_time(clock(1), clock(2))
        line = 'Starting on '//trim(adjustl(clock(1) (1:8)))//' at '//trim(adjustl(clock(2)))
        call TLab_Write_ASCII(lfile, line)

        line = 'Git-hash '//GITHASH
        call TLab_Write_ASCII(lfile, line)

        line = 'Git-branch '//GITBRANCH
        call TLab_Write_ASCII(lfile, line)

        !#####################################################################
        ! Inititalize OpenMP mode
#ifdef USE_OPENMP
        TLab_OMP_numThreads = omp_get_max_threads()
        write (line, *) TLab_OMP_numThreads
        line = 'Number of OMP threads '//trim(adjustl(line))
        call TLab_Write_ASCII(lfile, line)

#else
        TLab_OMP_numThreads = 1

#endif

        return
    end subroutine TLab_Start

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Stop(error_code)
        integer, intent(in) :: error_code

        ! ###################################################################
        if (error_code /= 0) then
            write (line, *) error_code
            line = 'Error code '//trim(adjustl(line))//'.'
            call TLab_Write_ASCII(efile, line)
        end if

        call GETARG(0, line)
        write (line, *) 'Finalizing program '//trim(adjustl(line))
        if (error_code == 0) then
            line = trim(adjustl(line))//' normally.'
        else
            line = trim(adjustl(line))//' abnormally. Check '//trim(adjustl(efile))
        end if
        call TLab_Write_ASCII(lfile, line)

#ifdef USE_MPI
        ims_time_max = MPI_WTIME()
        write (line, fmt_r) ims_time_max - ims_time_min
        line = 'Time elapse ....................: '//trim(adjustl(line))
        call TLab_Write_ASCII(lfile, line)

#ifdef PROFILE_ON
        write (line, fmt_r) ims_time_trans
        line = 'Time in array transposition ....: '//trim(ADJUST(line))
        call TLab_Write_ASCII(lfile, line)
#endif

#endif

        call TLab_Write_ASCII(lfile, '########################################')

#ifdef USE_MPI
        if (error_code == 0) then
            call MPI_FINALIZE(ims_err)
        else
            call MPI_Abort(MPI_COMM_WORLD, error_code, ims_err)
        end if
#else
        stop
#endif

        return
    end subroutine TLab_Stop

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Write_ASCII(file, locLine, flag_all)
        character*(*), intent(in) :: file, locLine
        logical, intent(in), optional :: flag_all

        ! -----------------------------------------------------------------------
        character*10 clock(2)

        ! #######################################################################
#ifdef USE_MPI
        if (ims_pro == 0 .or. present(flag_all)) then
#endif

            if (imode_verbosity > 0) then

                open (UNIT=22, FILE=file, STATUS='unknown', POSITION='APPEND')
                if (imode_verbosity == 1) then
                    write (22, '(a)') trim(adjustl(locLine))
                else if (imode_verbosity == 2) then
                    call date_and_time(clock(1), clock(2))
                    write (22, '(a)') '['//trim(adjustr(clock(2)))//'] '//trim(adjustl(locLine))
                end if
                close (22)

            end if

#ifndef PARALLEL
            if (file == efile) then
                write (*, *) trim(adjustl(locLine))
            end if
#endif

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine TLab_Write_ASCII

end module TLab_WorkFlow
