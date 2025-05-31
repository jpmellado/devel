#include "tlab_error.h"

module Microphysics
    use TLab_Constants, only: wp, wi, pi_wp, efile, MAX_VARS, MAX_PARS
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC, settling
    use TLab_Memory, only: inb_scal_array
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_Anelastic, only: dsmooth
    use Thermo_AirWater, only: inb_scal_ql
    use OPR_Partial, only: OPR_Partial_Z, OPR_P1
    ! use OPR_ODES
    implicit none
    private

    public :: Microphysics_Initialize
    public :: Microphysics_Sedimentation

    ! -------------------------------------------------------------------
    type microphysics_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical :: active(MAX_VARS) = .false.   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type microphysics_dt

    type(microphysics_dt), public, protected :: sedimentationProps          ! Microphysics parameters
    integer, parameter :: TYPE_SED_NONE = 0
    integer, parameter :: TYPE_SED_AIRWATER = 1
    integer, parameter :: TYPE_SED_AIRWATERSIMPLIFIED = 2

    type(microphysics_dt), public, protected :: evaporationProps            ! Microphysics parameters
    integer, parameter :: TYPE_EVA_NONE = 0
    integer, parameter, public :: TYPE_EVA_EQUILIBRIUM = 1

contains
!########################################################################
!########################################################################
    subroutine Microphysics_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer(wi) idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        ! ###################################################################
        block = 'Evaporation'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; evaporationProps%type = TYPE_EVA_NONE
        elseif (trim(adjustl(sRes)) == 'equilibrium') then; evaporationProps%type = TYPE_EVA_EQUILIBRIUM
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (evaporationProps%type /= TYPE_SED_NONE) then
            evaporationProps%parameters(:) = 1.0_wp         ! default values
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = MAX_PARS
                call LIST_REAL(sRes, idummy, evaporationProps%parameters)
            end if

        end if

        ! -------------------------------------------------------------------
        ! Initialize
        if (imixture == MIXT_TYPE_AIRWATER .and. evaporationProps%type == TYPE_EVA_EQUILIBRIUM) then
            inb_scal_array = inb_scal_array + 1         ! space for ql as diagnostic variable.

!            use Thermo_Base, only: dsmooth, NEWTONRAPHSON_ERROR
            dsmooth = evaporationProps%parameters(1)    ! Smooth factor for thermodynamic partial derivatives

        end if

        ! ###################################################################
        ! Read
        block = 'Sedimentation'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')
        call TLab_Write_ASCII(bakfile, '#Exponent=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; sedimentationProps%type = TYPE_SED_NONE
        elseif (trim(adjustl(sRes)) == 'airwater') then; sedimentationProps%type = TYPE_SED_AIRWATER
        elseif (trim(adjustl(sRes)) == 'airwatersimplified') then; sedimentationProps%type = TYPE_SED_AIRWATERSIMPLIFIED
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (sedimentationProps%type /= TYPE_SED_NONE) then
            sedimentationProps%parameters(:) = 1.0_wp           ! default values
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = MAX_PARS
                call LIST_REAL(sRes, idummy, sedimentationProps%parameters)
            end if

            if (imixture == MIXT_TYPE_AIRWATER) then
                sedimentationProps%active(:) = .true.           ! All scalars are affected
                call ScanFile_Real(bakfile, inifile, block, 'Exponent', '0.0', sedimentationProps%auxiliar(1))
            end if

        end if

        ! -------------------------------------------------------------------
        ! Initialize
        sedimentationProps%scalar = inb_scal_array      ! By default, sedimentation is caused by last scalar
        if (imixture == MIXT_TYPE_AIRWATER) then
            sedimentationProps%scalar = inb_scal_ql     ! sedimentation is caused by liquid
        end if

        if (sedimentationProps%type /= TYPE_SED_NONE) then
            if (settling > 0.0_wp) then
                sedimentationProps%parameters = sedimentationProps%parameters*settling ! adding the settling number in the parameter definitions
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Settling number must be nonzero if sedimentation is retained.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        return
    end subroutine Microphysics_Initialize

!########################################################################
!########################################################################
    subroutine Microphysics_Sedimentation(locProps, nx, ny, nz, is, g, s, source, tmp1, flux)
        use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_OutPlace, Thermo_Anelastic_StaticL
        use FDM, only: fdm_dt
        type(microphysics_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        target s, source

        ! -----------------------------------------------------------------------
        real(wp) dummy, exponent

        real(wp), pointer :: s_active(:) => null()

        !########################################################################
        exponent = locProps%auxiliar(1)
        dummy = 1.0_wp + exponent

        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, rbackground, s(:, locProps%scalar(is)), source)
            s_active => source
        else
            s_active => s(:, locProps%scalar(is))
        end if

        select case (locProps%type)
        case (TYPE_SED_AIRWATER)
            select case (is)
            case (2, 3)         ! q_t, q_l
                if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                    tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*(s_active**dummy)
                else
                    tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*s_active
                end if

            case default        ! energy variables
                call Thermo_Anelastic_StaticL(nx, ny, nz, s, tmp1)
                if (exponent > 0.0_wp) then
                    tmp1 = locProps%parameters(is)*tmp1*(s_active**dummy)
                else
                    tmp1 = locProps%parameters(is)*tmp1*s_active
                end if

            end select

            call OPR_Partial_Z(OPR_P1, nx, ny, nz, g, tmp1, source)
            if (present(flux)) flux = -tmp1

        case (TYPE_SED_AIRWATERSIMPLIFIED)
            ! if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
            !     tmp1 = locProps%parameters(is)*(s_active**dummy)
            ! else
            !     tmp1 = locProps%parameters(is)*s_active
            ! end if

            ! call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g, tmp1, source)
            ! if (present(flux)) flux = -tmp1

            ! the previous formulation yields oscillations at sharp gradients
            call OPR_Partial_Z(OPR_P1, nx, ny, nz, g, s_active, tmp1)
            if (exponent > 0.0_wp) tmp1 = tmp1*(s_active**exponent)
            source = locProps%parameters(is)*dummy*tmp1

            if (present(flux)) flux = -locProps%parameters(is)*(s_active**dummy)

        end select

! ###################################################################
        nullify (s_active)

    end subroutine Microphysics_Sedimentation

! !########################################################################
! !########################################################################
!     subroutine Microphysics_Evaporation()

!         return
!     end subroutine Microphysics_Evaporation

end module
