#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Reading general data from file tlab.ini, setting up general parameters
!########################################################################
subroutine TLab_Initialize_Parameters(inifile)
    use TLab_Constants, only: wp, wi, lfile, efile, wfile, MajorVersion, MinorVersion
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_WorkFlow, only: imode_sim, imode_verbosity
    use TLab_WorkFlow, only: flow_on, scal_on, stagger_on
    use TLab_Memory, only: imax, jmax, kmax
    use IO_Fields, only: io_fileformat, io_datatype, IO_MPIIO, IO_NETCDF, IO_NOFILE, IO_TYPE_DOUBLE, IO_TYPE_SINGLE
    implicit none

    character(len=*), intent(in) :: inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*32 bakfile
    integer(wi) idummy

! ###################################################################
    bakfile = trim(adjustl(inifile))//'.bak'

    call TLab_Write_ASCII(lfile, 'Reading global input data.')

! ###################################################################
! Version Checking
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#[Version]')
    call TLab_Write_ASCII(bakfile, '#Major=<mayor version number>')
    call TLab_Write_ASCII(bakfile, '#Minor=<minor version number>')

    call ScanFile_Int(bakfile, inifile, 'Version', 'Major', '0', idummy)
    if (MajorVersion /= idummy) then
        call TLab_Write_ASCII(efile, __FILE__//'. Major version error.')
        call TLab_Stop(DNS_ERROR_VERSION)
    end if
    call ScanFile_Int(bakfile, inifile, 'Version', 'Minor', '0', idummy)
    if (MinorVersion /= idummy) then
        write (sRes, '(I5)') MinorVersion
        call TLab_Write_ASCII(wfile, 'DNS_REAL_GLOBAL. Minor version warning. Expected : '//sRes)
    end if

! ###################################################################
! Global information
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Main]')
    call TLab_Write_ASCII(bakfile, '#FileFormat=<mpiio/NetCDF/None>')
    call TLab_Write_ASCII(bakfile, '#FileType=<Double/Single>')
    call TLab_Write_ASCII(bakfile, '#VerbosityLevel=<0/1/2>')
    call TLab_Write_ASCII(bakfile, '#CalculateFlow=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#CalculateScalar=<yes/no>')

    ! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, inifile, 'Main', 'FileFormat', 'MpiIO', sRes)
    if (trim(adjustl(sRes)) == 'mpiio') then; io_fileformat = IO_MPIIO
    elseif (trim(adjustl(sRes)) == 'netcdf') then; io_fileformat = IO_NETCDF
    elseif (trim(adjustl(sRes)) == 'none') then; io_fileformat = IO_NOFILE
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong Main.FileFormat.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'FileType', 'Double', sRes)
    if (trim(adjustl(sRes)) == 'double') then; io_datatype = IO_TYPE_DOUBLE
    elseif (trim(adjustl(sRes)) == 'single') then; io_datatype = IO_TYPE_SINGLE
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong Main.FileType.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end if

    call ScanFile_Int(bakfile, inifile, 'Main', 'VerbosityLevel', '1', imode_verbosity)

    call ScanFile_Char(bakfile, inifile, 'Main', 'CalculateFlow', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; flow_on = .true.
    elseif (trim(adjustl(sRes)) == 'no') then; flow_on = .false.
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Entry Main.CalculateFlow must be yes or no')
        call TLab_Stop(DNS_ERROR_CALCFLOW)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'CalculateScalar', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; scal_on = .true.
    elseif (trim(adjustl(sRes)) == 'no') then; scal_on = .false.
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Entry Main.CalculateScalar must be yes or no')
        call TLab_Stop(DNS_ERROR_CALCSCALAR)
    end if

! ###################################################################
! Pressure staggering
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Staggering]')
    call TLab_Write_ASCII(bakfile, '#StaggerHorizontalPressure=<yes/no>')

    call ScanFile_Char(bakfile, inifile, 'Staggering', 'StaggerHorizontalPressure', 'no', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; stagger_on = .true.; call TLab_Write_ASCII(lfile, 'Horizontal staggering of the pressure along Ox and Oz.')
    elseif (trim(adjustl(sRes)) == 'no') then; stagger_on = .false.
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Entry Main. StaggerHorizontalPressure must be yes or no')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

! ###################################################################
! Grid size
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Grid]')
    call TLab_Write_ASCII(bakfile, '#Imax=<imax>')
    call TLab_Write_ASCII(bakfile, '#Jmax=<jmax>')
    call TLab_Write_ASCII(bakfile, '#Kmax=<kmax>')

    call ScanFile_Int(bakfile, inifile, 'Grid', 'Imax', '0', imax)
    call ScanFile_Int(bakfile, inifile, 'Grid', 'Jmax', '0', jmax)
    call ScanFile_Int(bakfile, inifile, 'Grid', 'Kmax', '0', kmax)

    return
end subroutine TLab_Initialize_Parameters
