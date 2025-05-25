module Thermo_Compressible
    use TLab_Constants, only: wp, wi
    use Thermo_Base, only: gama0, nondimensional
    use Thermo_Base, only: THERMO_R, THERMO_PSAT
    use Thermo_AirWater, only: PREF_1000
    ! use NavierStokes, only: mach
    implicit none
    private

    public :: Thermo_Compressible_Initialize

    ! Nondimensional formulation
    !                                       A dimensional formulation can be imposed by setting RRATIO=CRATIO_INV=1
    real(wp), public :: RRATIO = 1.0_wp     ! 1/(gama0 mach^2) = R0/(U0^2/T0)
    real(wp), public :: CRATIO_INV = 1.0_wp ! (gamma0-1)*mach^2 = (U0^2/T0)/Cp0
    real(wp), public :: RRATIO_INV = 1.0_wp ! gama0 mach^2 = (U0^2/T0)/R0, inverse of RRATIO to save computational time in some routines

contains
    !########################################################################
    !########################################################################
    subroutine Thermo_Compressible_Initialize()

        ! Parameters in the evolution equations
        ! Compressible formulations use p nondimensionalized by reference dynamic pressure rho_0 U_0^2
        if (nondimensional) then
            ! RRATIO = 1.0_wp/(gama0*mach*mach)       ! (R_0T_0)/U_0^2 = p_0/(rho_0U_0^2), a scaled reference pressure
            ! CRATIO_INV = (gama0 - 1.0_wp)*mach*mach

            PREF_1000 = PREF_1000*RRATIO            ! scaling by dynamic reference pressure U0^2/T0 for compressible mode
            THERMO_PSAT(:) = THERMO_PSAT(:)*RRATIO
            THERMO_R(:) = THERMO_R(:)*RRATIO

            RRATIO_INV = 1.0_wp/RRATIO

        ! else
        !     if (imixture == MIXT_TYPE_NONE) then
        !         call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Single species formulation must be nondimensional.')
        !         call TLab_Stop(DNS_ERROR_OPTION)
        !     end if

        end if

        return
    end subroutine Thermo_Compressible_Initialize

end module Thermo_Compressible
