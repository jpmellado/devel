module Thermo_AirWater
    use TLab_Constants, only: wi, wp
    use Thermo_Base, only: imixture
    use Thermo_Base, only: MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER
    use Thermo_Base, only: THERMO_R, THERMO_CP
    use Thermo_Base, only: PREF, nondimensional
    use TLab_Memory, only: inb_scal_array
    implicit none
    private

    public :: Thermo_AirWater_Initialize
    real(wp), public :: Rv, Rd, Rdv, Cd, Cl, Cdv, Cvl, Cdl, Lv0, Ld, Lv, Ldv, Lvl, Ldl, rd_ov_rv

    ! Reference to calculate potential temperatures
    real(wp), public :: PREF_1000 = 1e5_wp          ! 1000 hPa

contains
    ! ###################################################################
    ! ###################################################################
    subroutine Thermo_AirWater_Initialize()

        ! -------------------------------------------------------------------
        ! diagnostic variables
        if (imixture == MIXT_TYPE_AIRWATER) then    ! assumes phase equilibrium
            inb_scal_array = inb_scal_array + 1     ! liquid water specific humidity, q_l
            inb_scal_array = inb_scal_array + 1     ! temperature, T
        end if

        ! -------------------------------------------------------------------
        ! Definitions for clarity in the airwater thermodynamics code
        Rv = THERMO_R(1)
        Rd = THERMO_R(2)
        Rdv = THERMO_R(1) - THERMO_R(2)
        rd_ov_rv = Rd/Rv

        Cd = THERMO_CP(1, 1, 2)
        Cl = THERMO_CP(1, 1, 3)
        Cdv = THERMO_CP(1, 1, 1) - THERMO_CP(1, 1, 2)
        Cvl = THERMO_CP(1, 1, 3) - THERMO_CP(1, 1, 1)
        Cdl = THERMO_CP(1, 1, 3) - THERMO_CP(1, 1, 2)
        Lv0 = -THERMO_CP(6, 1, 3)
        Ld = THERMO_CP(6, 1, 2)
        Lv = THERMO_CP(6, 1, 1)
        Ldv = THERMO_CP(6, 1, 1) - THERMO_CP(6, 1, 2)
        Lvl = THERMO_CP(6, 1, 3) - THERMO_CP(6, 1, 1)
        Ldl = THERMO_CP(6, 1, 3) - THERMO_CP(6, 1, 2)

        return
    end subroutine Thermo_AirWater_Initialize

end module Thermo_AirWater
