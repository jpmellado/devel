#include "tlab_error.h"

!# Mixture:
!# NSP               # of species in the mixture (NSP>=inb_scal)
!# inb_scal          # of prognostic scalars (transported during the simulation and saved)
!# inb_scal_array    # of scalars in array s (prognositc+diagnostic, normally = inb_scal)
!#
!# General multispecies formulation retains NSP-1 species in scalar array s, Y_i=s_i, the last one is obtained by sum Y_i=1
!# There might be additional scalars at the end of array s, e.g., conserved scalars or diagnostic variables
!# Multispecies admits the case Y_i=f_i(s_j), s_j could be conserved scalar, e.g. using total water or mixture fraction
!#
!# Saturation pressure assumes that reference rho_0 in non-dimensionalization is such that reference pressure is 1 bar.

module Thermodynamics
    use TLab_Constants, only: wi, wp
    use TLab_Constants, only: efile, lfile, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thermo_Base
    use Thermo_AirWater
    use Thermo_Anelastic
    use Thermo_Compressible
    implicit none
    private

    public :: Thermo_Initialize

    integer, public :: imode_thermo
    integer, parameter, public :: THERMO_TYPE_NONE = 0
    integer, parameter, public :: THERMO_TYPE_ANELASTIC = 1
    integer, parameter, public :: THERMO_TYPE_COMPRESSIBLE = 2

contains
    !########################################################################
    !########################################################################
    subroutine Thermo_Initialize(inifile)
        character(len=*), intent(in), optional :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr, str
        character(len=512) sRes

        integer(wi) icp, is, im, ipsat
        real(wp) CPREF, RREF

        !########################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Thermodynamics'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<compressible/anelastic>')
        call TLab_Write_ASCII(bakfile, '#Gamma=<value>')
        call TLab_Write_ASCII(bakfile, '#Mixture=<value>')
        call TLab_Write_ASCII(bakfile, '#Nondimensional=<yes,no>')
        ! call TLab_Write_ASCII(bakfile, '#SmoothFactor=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; imode_thermo = THERMO_TYPE_NONE
        else if (trim(adjustl(sRes)) == 'anelastic') then; imode_thermo = THERMO_TYPE_ANELASTIC
        else if (trim(adjustl(sRes)) == 'compressible') then; imode_thermo = THERMO_TYPE_COMPRESSIBLE
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'Mixture', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; imixture = MIXT_TYPE_NONE
        else if (trim(adjustl(sRes)) == 'air') then; imixture = MIXT_TYPE_AIR
        else if (trim(adjustl(sRes)) == 'airvapor') then; imixture = MIXT_TYPE_AIRVAPOR
        else if (trim(adjustl(sRes)) == 'airwater') then; imixture = MIXT_TYPE_AIRWATER
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Mixture.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'Nondimensional', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') then; nondimensional = .true.
        else if (trim(adjustl(sRes)) == 'no') then; nondimensional = .false.
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Nondimensional.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Real(bakfile, inifile, block, 'HeatCapacityRatio', '1.4', gamma0)

        ! call ScanFile_Real(bakfile, inifile, block, 'SmoothFactor', '0.1', dsmooth)

        !########################################################################
        ! Initialize
        call Thermo_Mixture_Initialize()

        RREF = THERMO_R(ISPREF)                 ! Reference value R_0

        CPREF = 0.0_wp                          ! Reference heat capacity Cp0, J /kg /K
        do icp = NCP, 1, -1
            CPREF = CPREF*TREF + THERMO_CP(icp, 2, ISPREF)
        end do

        if (imixture /= MIXT_TYPE_NONE) then    ! Reference heat capacity ratio; otherwise, gamma0 is read in tlab.ini
            gamma0 = CPREF/(CPREF - RREF)
        end if

        ! -------------------------------------------------------------------
        ! Nondimensionalization
        if (nondimensional) then
            ! Thermal equation of state
            THERMO_R(1:NSP) = THERMO_R(1:NSP)/RREF      ! Normalized gas constants (Inverse of molar masses)

            ! Caloric equations of state
            do is = 1, NSP
                do im = 1, 2
                    THERMO_CP(6, im, is) = THERMO_CP(6, im, is)/(CPREF*TREF)    ! Formation enthalpy
                    THERMO_CP(7, im, is) = THERMO_CP(7, im, is)/CPREF           ! Formation entropy
                    do icp = 1, NCP
                        THERMO_CP(icp, im, is) = THERMO_CP(icp, im, is)*(TREF**(icp - 1))/CPREF
                    end do
                end do
            end do
            THERMO_TLIM = THERMO_TLIM/TREF              ! Temperature limits for polynomial fits to cp

            PREF_1000 = PREF_1000/PREF                  ! 1000 hPa, reference to calculate potential temperatures

            THERMO_PSAT(:) = THERMO_PSAT(:)/PREF        ! Saturation vapor pressure
            do ipsat = 1, NPSAT
                THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)*(TREF**(ipsat - 1))
            end do

        end if

        ! -------------------------------------------------------------------
        select case (imode_thermo)
        case (THERMO_TYPE_COMPRESSIBLE)
            call Thermo_Compressible_Initialize()

        case (THERMO_TYPE_ANELASTIC)
            call Thermo_Anelastic_Initialize(inifile)

        end select

        if (any([MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER] == imixture)) &
            call Thermo_AirWater_Initialize()

        ! ###################################################################
        ! Output
        ! ###################################################################
        if (imixture /= MIXT_TYPE_NONE) then

            call TLab_Write_ASCII(lfile, 'Thermodynamic properties have been initialized.')
            do is = 1, NSP
                write (str, *) is; str = 'Setting Species'//trim(adjustl(str))//'='//trim(adjustl(THERMO_SPNAME(is)))
                call TLab_Write_ASCII(lfile, str)
            end do

            write (str, fmt_r) TREF
            call TLab_Write_ASCII(lfile, 'Setting TRef = '//trim(adjustl(str)))

            write (str, fmt_r) PREF
            call TLab_Write_ASCII(lfile, 'Setting pRef = '//trim(adjustl(str)))

            write (str, fmt_r) PREF/(RREF*TREF)
            call TLab_Write_ASCII(lfile, 'Setting rhoRef = '//trim(adjustl(str)))

            write (str, fmt_r) RREF
            call TLab_Write_ASCII(lfile, 'Setting RRef = '//trim(adjustl(str)))

            write (str, fmt_r) CPREF
            call TLab_Write_ASCII(lfile, 'Setting CpRef = '//trim(adjustl(str)))

            write (str, fmt_r) gamma0
            call TLab_Write_ASCII(lfile, 'Setting Gamma0 = '//trim(adjustl(str)))

        end if

        return
    end subroutine Thermo_Initialize

end module Thermodynamics
