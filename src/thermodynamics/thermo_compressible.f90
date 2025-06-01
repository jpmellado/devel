module Thermo_Compressible
    use TLab_Constants, only: wp, wi
    use Thermo_Base, only: gamma0, nondimensional
    use Thermo_Base, only: THERMO_R, THERMO_PSAT
    use Thermo_AirWater, only: PREF_1000
    ! use NavierStokes, only: mach
    implicit none
    private

    public :: Thermo_Compressible_Initialize

    ! Nondimensional formulation
    !                                       A dimensional formulation can be imposed by setting RRATIO=CRATIO_INV=1
    real(wp), public :: RRATIO = 1.0_wp     ! 1/(gamma0 mach^2) = R0/(U0^2/T0)
    real(wp), public :: CRATIO_INV = 1.0_wp ! (gamma0-1)*mach^2 = (U0^2/T0)/Cp0
    real(wp), public :: RRATIO_INV = 1.0_wp ! gamma0 mach^2 = (U0^2/T0)/R0, inverse of RRATIO to save computational time in some routines

    ! pointers
    real(wp), pointer, public :: e(:, :, :) => null()
    real(wp), pointer, public :: rho(:, :, :) => null()
    real(wp), pointer, public :: p(:, :, :) => null()
    real(wp), pointer, public :: T(:, :, :) => null()
    real(wp), pointer, public :: vis(:, :, :) => null()

    contains
    !########################################################################
    !########################################################################
    subroutine Thermo_Compressible_Initialize()

        ! Parameters in the evolution equations
        ! Compressible formulations use p nondimensionalized by reference dynamic pressure rho_0 U_0^2
        if (nondimensional) then
            ! RRATIO = 1.0_wp/(gamma0*mach*mach)       ! (R_0T_0)/U_0^2 = p_0/(rho_0U_0^2), a scaled reference pressure
            ! CRATIO_INV = (gamma0 - 1.0_wp)*mach*mach

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

    ! !########################################################################
    ! !########################################################################
    ! subroutine Thermo_Compressible_Memory()
    !     use TLab_Memory, only: isize_field, imax, jmax, kmax
    !     use TLab_Arrays, only: q

    !     integer(wi) idummy(2)

    !     idummy = shape(q)
    !     ! Pointers
    !     if (idummy(2) >= 4) e(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 4)
    !     if (idummy(2) >= 5) rho(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 5)
    !     if (idummy(2) >= 6) p(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 6)
    !     if (idummy(2) >= 7) T(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 7)
    !     if (idummy(2) >= 8) vis(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 8)

    !     return
    ! end subroutine Thermo_Compressible_Memory

end module Thermo_Compressible
