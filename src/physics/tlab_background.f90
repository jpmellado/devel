#include "tlab_error.h"

! information to set up bcs, ics, and reference background profiles

module Tlab_Background
    use TLab_Constants, only: wp, wi, efile, lfile, wfile, MAX_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Profiles, only: profiles_dt, Profiles_Calculate
    implicit none
    private

    public :: TLab_Initialize_Background

    type(profiles_dt), public :: qbg(3)                     ! Velocity
    type(profiles_dt), public :: sbg(MAX_VARS)              ! Scalars
    type(profiles_dt), public :: pbg, rbg, tbg, hbg         ! Pressure, density, temperature, enthalpy

    ! background, reference profiles
    real(wp), allocatable :: sbackground(:, :)              ! Scalars

contains
!########################################################################
!# Initialize data of reference profiles
!########################################################################
    subroutine TLab_Initialize_Background(inifile)
        use TLab_Arrays, only: wrk1d
        use TLab_Memory, only: inb_scal, inb_scal_array
        use TLab_Grid, only: z
        use FDM, only: fdm_Int0
        use NavierStokes, only: schmidt
        use Thermodynamics, only: imode_thermo, THERMO_TYPE_ANELASTIC
        use Thermo_Anelastic
        use Profiles, only: Profiles_ReadBlock
        use Profiles, only: PROFILE_NONE, PROFILE_EKMAN_U, PROFILE_EKMAN_U_P, PROFILE_EKMAN_V
        use Gravity

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile, block
        ! character(len=512) sRes
        character(len=64) lstr

        ! real(wp) ploc(2), rloc(2), Tloc(2), sloc(2, 1:MAX_VARS)

        integer(wi) is, k!, idummy

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        ! -----------------------------------------------------------------------
        block = 'Scalar'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        do is = 1, MAX_VARS
            write (lstr, *) is
            call Profiles_ReadBlock(bakfile, inifile, block, 'Scalar'//trim(adjustl(lstr)), sbg(is))
        end do

        ! -----------------------------------------------------------------------
        block = 'Flow'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityX', qbg(1))
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityY', qbg(2))
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityZ', qbg(3))

        ! Consistency check
        if (any([PROFILE_EKMAN_U, PROFILE_EKMAN_U_P] == qbg(1)%type)) then
            qbg(3) = qbg(1)
            qbg(3)%type = PROFILE_EKMAN_V
        end if

        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Pressure', pbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Density', rbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Temperature', tbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Enthalpy', hbg)

        ! ! consistency check; two and only two are given TO BE CHECKED BECAUSE PROFILE_NONE is used as constant profile
        ! if (imode_thermo == THERMO_TYPE_COMPRESSIBLE) then
        !     idummy = 0
        !     if (pbg%type == PROFILE_NONE) idummy = idummy + 1
        !     if (rbg%type == PROFILE_NONE) idummy = idummy + 1
        !     if (tbg%type == PROFILE_NONE) idummy = idummy + 1
        !     if (hbg%type == PROFILE_NONE) idummy = idummy + 1
        !     if (idummy /= 2) then
        !         call TLab_Write_ASCII(efile, __FILE__//'. Specify only 2 thermodynamic profiles.')
        !         call TLab_Stop(DNS_ERROR_OPTION)
        !     end if
        ! end if

        ! ! #######################################################################
        ! ! mean_rho and delta_rho need to be defined, because of old version.
        ! if (imode_thermo == THERMO_TYPE_COMPRESSIBLE) then
        !     if (rbg%type == PROFILE_NONE .and. tbg%type /= PROFILE_NONE) then
        !         rbg = tbg

        !         Tloc(1) = tbg%mean + 0.5_wp*tbg%delta
        !         Tloc(2) = tbg%mean - 0.5_wp*tbg%delta
        !         ploc(1:2) = pbg%mean
        !         sloc(1, :) = sbg(:)%mean + 0.5_wp*sbg(:)%delta
        !         sloc(2, :) = sbg(:)%mean - 0.5_wp*sbg(:)%delta
        !         call THERMO_THERMAL_DENSITY(2, sloc, ploc, Tloc, rloc)
        !         rbg%mean = 0.5_wp*(rloc(1) + rloc(2))
        !         rbg%delta = rloc(1) - rloc(2)

        !     else if (rbg%type == PROFILE_NONE .and. hbg%type /= PROFILE_NONE) then
        !         rbg%mean = 1.0_wp ! to be done; I need a nonzero value for dns_control.

        !     else
        !         tbg = rbg

        !         rloc(1) = rbg%mean + 0.5_wp*rbg%delta
        !         rloc(2) = rbg%mean - 0.5_wp*rbg%delta
        !         ploc(1:2) = pbg%mean
        !         sloc(1, :) = sbg(:)%mean + 0.5_wp*sbg(:)%delta
        !         sloc(2, :) = sbg(:)%mean - 0.5_wp*sbg(:)%delta
        !         call THERMO_THERMAL_TEMPERATURE(2, sloc, ploc, rloc, Tloc)
        !         tbg%mean = 0.5_wp*(Tloc(1) + Tloc(2))
        !         tbg%delta = Tloc(1) - Tloc(2)

        !     end if
        ! end if

        do is = 1, size(qbg)
            if (qbg(is)%relative) qbg(is)%zmean = z%nodes(1) + z%scale*qbg(is)%zmean_rel
        end do
        if (pbg%relative) pbg%zmean = z%nodes(1) + z%scale*pbg%zmean_rel
        if (rbg%relative) rbg%zmean = z%nodes(1) + z%scale*rbg%zmean_rel
        if (tbg%relative) tbg%zmean = z%nodes(1) + z%scale*tbg%zmean_rel
        if (hbg%relative) hbg%zmean = z%nodes(1) + z%scale*hbg%zmean_rel
        do is = 1, size(sbg)
            if (sbg(is)%relative) sbg(is)%zmean = z%nodes(1) + z%scale*sbg(is)%zmean_rel
        end do

        ! #######################################################################
        ! -----------------------------------------------------------------------
        ! Construct reference  profiles
        allocate (sbackground(z%size, inb_scal_array))   ! scalar profiles

        do is = 1, inb_scal
            do k = 1, z%size
                sbackground(k, is) = Profiles_Calculate(sbg(is), z%nodes(k))
            end do
        end do

        if (imode_thermo == THERMO_TYPE_ANELASTIC) then     ! thermodynamic profiles
            allocate (epbackground(z%size))
            allocate (tbackground(z%size))
            allocate (pbackground(z%size))
            allocate (rbackground(z%size))
            allocate (ribackground(z%size))

            call Gravity_Hydrostatic_Enthalpy(fdm_Int0, z%nodes(:), sbackground, epbackground, tbackground, pbackground, pbg%zmean, pbg%mean, wrk1d)

            call Thermo_Anelastic_Rho(1, 1, z%size, sbackground, rbackground, wrk1d)
            ribackground = 1.0_wp/rbackground

        end if

        if (any(gravityProps%active(1:3))) then
            allocate (bbackground(z%size))                   ! buoyancy profiles
            bbackground(:) = 0.0_wp

            if (gravityProps%active(3)) then
                call Gravity_Source(gravityProps, 1, 1, z%size, sbackground(:, 1), wrk1d)
                bbackground(1:z%size) = wrk1d(1:z%size, 1)
            end if

        end if

        ! -----------------------------------------------------------------------
        ! Add diagnostic fields to reference profile data, if any
        do is = inb_scal + 1, inb_scal_array ! Add diagnostic fields, if any
            sbg(is) = sbg(1)
            schmidt(is) = schmidt(1)
        end do

        return
    end subroutine TLab_Initialize_Background

end module Tlab_Background
