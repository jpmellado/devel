#include "tlab_error.h"

module Buffer
    use TLab_Constants, only: wp, wi, MAX_PARS
    use TLab_Constants, only: efile, lfile
    use TLab_Constants, only: tag_flow, tag_scal, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: inb_flow, inb_scal
    implicit none
    private

    integer, parameter, public :: BUFFER_TYPE_NONE = 0
    integer, parameter, public :: BUFFER_TYPE_NUDGE = 1
    ! integer, parameter :: BUFFER_TYPE_FILTER = 2
    integer(wi), public :: bufferType = BUFFER_TYPE_NONE

    public :: Buffer_Initialize
    public :: Buffer_Nudge

    ! -------------------------------------------------------------------
    logical, public :: bufferLoad

    type :: buffer_dt
        sequence
        integer type
        integer(wi) size, offset
        real(wp), allocatable :: tau(:)         ! relaxation timescale
        real(wp), allocatable :: ref(:, :, :)   ! reference field
        integer(wi) form                        ! form of function of relaxation term
        real(wp) :: parameters(MAX_PARS)
    end type buffer_dt
    type(buffer_dt), allocatable :: bufferFlowKmin(:), bufferFlowKmax(:)
    type(buffer_dt), allocatable :: bufferScalKmin(:), bufferScalKmax(:)

    integer(wi), parameter :: FORM_POWER_MIN = 1
    integer(wi), parameter :: FORM_POWER_MAX = 2

contains
    !########################################################################
    !########################################################################
    subroutine Buffer_Initialize(inifile)
        use TLab_Pointers_3D, only: p_q, p_s
        use TLab_Grid, only: z

        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, name
        character(len=128) eStr
        character(len=512) sRes
        integer is, iq

        !########################################################################
        ! Read
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'BufferZone'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/relaxation/filter/both>')
        call TLab_Write_ASCII(bakfile, '#LoadBuffer=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#PointsKmin=<value>')
        call TLab_Write_ASCII(bakfile, '#PointsKmax=<value>')
        call TLab_Write_ASCII(bakfile, '#ParametersKmin=<values>')
        call TLab_Write_ASCII(bakfile, '#ParametersKmax=<values>')
        ! call TLab_Write_ASCII(bakfile, '#PointsImin=<value>')         ! so far, only in K direction
        ! call TLab_Write_ASCII(bakfile, '#PointsImax=<value>')
        ! call TLab_Write_ASCII(bakfile, '#ParametersImin=<values>')
        ! call TLab_Write_ASCII(bakfile, '#ParametersImax=<values>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; bufferType = BUFFER_TYPE_NONE
        else if (trim(adjustl(sRes)) == 'relaxation') then; bufferType = BUFFER_TYPE_NUDGE
            ! else if (trim(adjustl(sRes)) == 'filter') then; bufferType = BUFFER_TYPE_FILTER
            ! else if (trim(adjustl(sRes)) == 'both') then; bufferType = BUFFER_TYPE_BOTH
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong Type option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'LoadBuffer', 'no', sRes)
        if (trim(adjustl(sRes)) == 'yes') then
            bufferLoad = .true.
        else
            bufferLoad = .false.
        end if

        if (bufferType /= BUFFER_TYPE_NONE) then
            allocate (bufferFlowKmin(1:inb_flow), bufferFlowKmax(1:inb_flow))
            bufferFlowKmin(:)%type = bufferType; bufferFlowKmax(:)%type = bufferType            ! So far, all the same
            bufferFlowKmin(:)%form = FORM_POWER_MIN; bufferFlowKmax(:)%form = FORM_POWER_MAX
            do iq = 1, inb_flow
                call Buffer_ReadBlock(bakfile, inifile, block, 'UKmin', bufferFlowKmin(iq))
                bufferFlowKmin(iq)%offset = 0

                call Buffer_ReadBlock(bakfile, inifile, block, 'UKmax', bufferFlowKmax(iq))
                bufferFlowKmax(iq)%offset = z%size - bufferFlowKmax(iq)%size

            end do

            allocate (bufferScalKmin(1:inb_scal), bufferScalKmax(1:inb_scal))
            bufferScalKmin(:)%type = bufferType; bufferScalKmax(:)%type = bufferType; ! So far, all the same
            bufferScalKmin(:)%form = FORM_POWER_MIN; bufferScalKmax(:)%form = FORM_POWER_MAX
            do is = 1, inb_scal
                call Buffer_ReadBlock(bakfile, inifile, block, 'SKmin', bufferScalKmin(is))
                bufferScalKmin(is)%offset = 0

                call Buffer_ReadBlock(bakfile, inifile, block, 'SKmax', bufferScalKmax(is))
                bufferScalKmax(is)%offset = z%size - bufferScalKmax(is)%size

            end do

        end if

        !########################################################################
        ! Initialize data
        if (bufferType /= BUFFER_TYPE_NONE) then
            do iq = 1, inb_flow
                write (name, *) iq; name = trim(adjustl(tag_flow))//'bcs.kmin.'//trim(adjustl(name))
                if (bufferFlowKmin(iq)%size > 0) call IniBlock_K(trim(adjustl(name)), bufferFlowKmin(iq), p_q(:, :, :, iq))
                write (name, *) iq; name = trim(adjustl(tag_flow))//'bcs.kmax.'//trim(adjustl(name))
                if (bufferFlowKmax(iq)%size > 0) call IniBlock_K(trim(adjustl(name)), bufferFlowKmax(iq), p_q(:, :, :, iq))
            end do

            do is = 1, inb_scal
                write (name, *) is; name = trim(adjustl(tag_scal))//'bcs.kmin.'//trim(adjustl(name))
                if (bufferScalKmin(is)%size > 0) call IniBlock_K(trim(adjustl(name)), bufferScalKmin(is), p_s(:, :, :, is))
                write (name, *) is; name = trim(adjustl(tag_scal))//'bcs.kmax.'//trim(adjustl(name))
                if (bufferScalKmax(is)%size > 0) call IniBlock_K(trim(adjustl(name)), bufferScalKmax(is), p_s(:, :, :, is))
            end do

        end if

        return
    end subroutine Buffer_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine Buffer_ReadBlock(bakfile, inifile, block, tag, locVar)
        character(len=*) bakfile, inifile, block, tag
        type(buffer_dt) locVar

        character(len=512) sRes
        character(len=128) eStr
        integer idummy

        ! ###################################################################
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call ScanFile_Int(bakfile, inifile, block, 'Points'//trim(adjustl(tag)), '0', locVar%size)

        if (locVar%size > 0) then
            call ScanFile_Char(bakfile, inifile, block, 'Parameters'//trim(adjustl(tag)), '1.0, 2.0', sRes)

            idummy = 2; call LIST_REAL(sRes, idummy, locVar%parameters)
            if (idummy /= 2) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong number of values in buffer parameters.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            if (locVar%parameters(1) == 0.0_wp) locVar%type = BUFFER_TYPE_NONE

        else
            locVar%type = BUFFER_TYPE_NONE

        end if

    end subroutine Buffer_ReadBlock

    ! ###################################################################
    ! ###################################################################
    subroutine IniBlock_K(tag, locVar, field)
        use Averages, only: AVG1V2D
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Grid, only: z
        use IO_Fields
        use TLab_Time, only: itime
#ifdef USE_MPI
        use mpi_f08, only: MPI_COMM_WORLD, MPI_REAL8
#endif
        character(len=*), intent(in) :: tag         ! File name information
        type(buffer_dt), intent(inout) :: locVar
        real(wp), intent(in) :: field(:, :, :)

        character(len=32) str, name
        character(len=128) line
        integer(wi) k, kglobal
        real(wp) dummy
        type(io_subarray_dt) :: io_subarrays
        integer(wi) io_sizes(5), idummy

        ! ###################################################################
        ! Reference fields
        allocate (locVar%ref(imax, jmax, locVar%size))

        io_subarrays%offset = 0
        io_subarrays%precision = wp
#ifdef USE_MPI
        io_subarrays%active = .true.
        io_subarrays%communicator = MPI_COMM_WORLD
        io_subarrays%subarray = IO_Create_Subarray_XOY(imax, jmax, locVar%size, MPI_REAL8)
#endif
        idummy = imax*jmax*locVar%size; io_sizes = (/idummy, 1, idummy, 1, 1/)

        name = trim(adjustl(tag))

        if (bufferLoad) then
            call IO_Read_Subarray(io_subarrays, name, [''], locVar%ref, io_sizes)

        else
            do k = 1, locVar%size
                kglobal = k + locVar%offset
                locVar%ref(:, :, k) = AVG1V2D(imax, jmax, kmax, kglobal, 1, field)
            end do

            idummy = index(name, '.', back=.true.)
            write (str, *) itime; name = trim(name(1:idummy - 1))//'.'//trim(adjustl(str))//trim(name(idummy:))
            call IO_Write_Subarray(io_subarrays, name, [''], locVar%ref, io_sizes)

        end if

        ! Control
        line = 'Checking bounds of field '//trim(adjustl(tag))//':'
        write (str, fmt_r) minval(locVar%ref(:, :, :)); line = trim(adjustl(line))//' '//trim(adjustl(str))//','
        write (str, fmt_r) maxval(locVar%ref(:, :, :)); line = trim(adjustl(line))//' '//trim(adjustl(str))//'.'
        call TLab_Write_ASCII(lfile, line)

        ! ###################################################################
        ! Inverse of relaxation time
        allocate (locVar%tau(locVar%size))

#define strength parameters(1)
#define sigma parameters(2)

        dummy = 1.0_wp/(z%nodes(locVar%offset + locVar%size) - z%nodes(locVar%offset + 1)) ! Inverse of segment length
        do k = 1, locVar%size
            kglobal = locVar%offset + k
            if (locVar%form == FORM_POWER_MAX) &
                locVar%tau(k) = locVar%strength*((z%nodes(kglobal) - z%nodes(locVar%offset + 1))*dummy)**locVar%sigma
            if (locVar%form == FORM_POWER_MIN) &
                locVar%tau(k) = locVar%strength*((z%nodes(locVar%offset + locVar%size) - z%nodes(kglobal))*dummy)**locVar%sigma
        end do

#undef strength
#undef sigma

        return
    end subroutine IniBlock_K

    ! ###################################################################
    ! ###################################################################
    subroutine Buffer_Nudge()
        use TLab_Pointers_3D, only: p_q, p_s
        use DNS_Arrays, only: p_hq, p_hs

        integer is, iq

        do iq = 1, inb_flow
            if (bufferFlowKmin(iq)%type == BUFFER_TYPE_NUDGE) call Nudge_K(bufferFlowKmin(iq), p_q(:, :, :, iq), p_hq(:, :, :, iq))
            if (bufferFlowKmax(iq)%type == BUFFER_TYPE_NUDGE) call Nudge_K(bufferFlowKmax(iq), p_q(:, :, :, iq), p_hq(:, :, :, iq))
        end do

        do is = 1, inb_scal
            if (bufferScalKmin(is)%type == BUFFER_TYPE_NUDGE) call Nudge_K(bufferScalKmin(is), p_s(:, :, :, is), p_hs(:, :, :, is))
            if (bufferScalKmax(is)%type == BUFFER_TYPE_NUDGE) call Nudge_K(bufferScalKmax(is), p_s(:, :, :, is), p_hs(:, :, :, is))
        end do

        return
    end subroutine Buffer_Nudge

    ! ###################################################################
    ! ###################################################################
    subroutine Nudge_K(locProps, s, hs)
        type(buffer_dt), intent(in) :: locProps
        real(wp), intent(in) :: s(:, :, :)
        real(wp), intent(out) :: hs(:, :, :)

        integer(wi) k, kglobal

        do k = 1, locProps%size
            kglobal = k + locProps%offset
            hs(:, :, kglobal) = hs(:, :, kglobal) - locProps%tau(k)*(s(:, :, kglobal) - locProps%ref(:, :, k))
        end do

        return
    end subroutine Nudge_K

end module Buffer
