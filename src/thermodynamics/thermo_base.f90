module Thermo_Base
    use TLab_Constants, only: wi, wp
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: Thermo_Mixture_Initialize
    public :: Thermo_Psat_Polynomial, Thermo_dPsat_Polynomial

    integer, public :: imixture
    integer, parameter, public :: MIXT_TYPE_NONE = 0
    integer, parameter, public :: MIXT_TYPE_AIR = 1
    integer, parameter, public :: MIXT_TYPE_AIRVAPOR = 2
    integer, parameter, public :: MIXT_TYPE_AIRWATER = 3
    integer, parameter, public :: MIXT_TYPE_CHEMKIN = -1

    integer, parameter, public :: MAX_NSP = 10          ! Maximum number of components (species) in a mixture
    integer(wi), public :: NSP = 0                      ! Number of components (species) in a mixture
    character(len=32), public :: THERMO_SPNAME(MAX_NSP) = ''

    ! Thermal data
    real(wp), public :: THERMO_R(MAX_NSP)               ! Gas constants

    ! Caloric data; Chemkin format
    integer(wi), parameter :: MAX_NCP = 7               ! Caloric data; cp(T), formation enthalpy, formation entropy
    integer(wi), public :: NCP                          ! Number of terms in polynomial fit to cp
    real(wp), public :: THERMO_CP(MAX_NCP, 2, MAX_NSP)  ! Polynomial coefficients. Last to terms are formation enthalpy and formation entropy
    !                                                   The second index indicates each of the 2 temperature intervals considered in the fit
    real(wp), public :: THERMO_TLIM(3, MAX_NSP)         ! Temperature limits of the two temperature intervals in the polynomial fits.

    ! Saturation vapor pressures
    integer(wi), parameter :: MAX_NPSAT = 10            ! Polynomial fit to saturation pressure
    integer(wi), public :: NPSAT
    real(wp), public :: THERMO_PSAT(MAX_NPSAT)

    ! ! Equilibrium calculations
    ! real(wp), public :: NEWTONRAPHSON_ERROR, dsmooth

    real(wp), public :: gamma0                           ! Specific heat ratio, Cp0/Cv0 = Cp0/(Cp0-R0)
    !                                                     For imixture=NONE, I only need gamma0 and can be set in tlab.ini
    !                                                     Otherwise, I need the thermodynamic data that is given in thermo_initialize, and gamma0 is derived.

    ! Nondimensional formulation
    logical, public :: nondimensional = .true.          ! consider nondimensional formulation
    ! Reference values for the nondimensionalization; in principle, these could be changed, if desired.
    integer, public :: ISPREF = 1                       ! Reference species for nondimensionalization
    real(wp), parameter, public :: TREF = 298.0_wp      ! K, reference temperature T_0
    real(wp), parameter, public :: PREF = 1e5_wp        ! Pa, reference pressure p_0
    !                                                     Reference density results from rho_0=p_0/(T_0R_0)

    ! -------------------------------------------------------------------
    real(wp), parameter :: RGAS = 8314_wp               ! J /kg /K, universal gas constant,

contains
    ! ###################################################################
    ! ###################################################################
    subroutine Thermo_Mixture_Initialize()
        real(wp) WGHT(MAX_NSP)                  ! Molar masses
        logical :: molar_data = .false.         ! Values are per unit mol or per unit mass
        integer is

        real(wp) TREF_LOC, HREF_LOC(MAX_NSP), SREF_LOC(MAX_NSP)

        real(wp) WRK1D_LOC(MAX_NPSAT), tmp1, tmp2
        integer(wi) ipsat, i, j

        ! ###################################################################
        ! Species tags
        ! Thermal equation, molar masses in kg/kmol
        ! ###################################################################
        select case (imixture)
            ! Iribarne and Godson, 1981
        case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2Ov'; WGHT(NSP) = 18.015_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Air'; WGHT(NSP) = 28.9644_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2Ol'; WGHT(NSP) = 18.015_wp
            ISPREF = 2

        end select

        THERMO_R(1:NSP) = RGAS/WGHT(1:NSP)              ! Specific gas constants, J /kg /K

        ! ###################################################################
        ! Caloric equations
        !
        ! General formulation is CHEMKIN format: 7-coefficient NASA polynomials (see Burcat&Ruscic):
        !
        ! THERMO_CP(i,im,k) = a_i of species k at
        !         im=1 high temperature
        !         im=2 low   temperature
        !
        ! C_{p,i} = \sum_1^5 a_i T^{i-1}
        ! h_i     = \sum_1^5 a_i T^i/i + a_6
        ! s_{T,i} = a_1 ln(T) + \sum_2^5 a_i T^{i-1}/(i-1) + a_7
        !
        ! i.e., dh_i = C_{p,i}dT and ds_{T,i}=C_{p,i}dT/T, where a_6 is
        ! related to the formation enthalpy, and a_7 to the formation entropy.
        ! HREF_LOC and SREF_LOC (at TREF_LOC) are used to fix last Cpi 6-7 coefficients.
        !
        ! Note that Burcat&Ruscic give values devided by R^0
        !
        ! NCP gives the number of coefficients a_i.
        ! The simplified situation NCP=1 assumes also only one range, which then gives C_p constant.
        ! This is used to expedite the calculation of T in the energy formulation.
        !
        ! The pressure contribution to the entropy still needs to be added
        !
        ! ###################################################################
        THERMO_CP(:, :, :) = 0.0_wp     ! Initialize to zero
        HREF_LOC(:) = 0.0_wp
        SREF_LOC(:) = 0.0_wp

        THERMO_TLIM(1, :) = 200.0_wp    ! Default T limits for the 2 intervals of polynomial fits
        THERMO_TLIM(2, :) = 5000.0_wp   ! These intervals are currently not used
        THERMO_TLIM(3, :) = 5000.0_wp

        NCP = 1                         ! Default order of heat capacity polynomial

        select case (imixture)
            ! Iribarne and Godson, 1981
        case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER)
            TREF_LOC = 273.15_wp

            molar_data = .false. ! this block gives data already in mass spefic values

            ! Enthalpy of Formation in J /kg
            HREF_LOC(1) = 1870.0_wp*TREF_LOC                 ! values s.t. THERMO_CP(6,im,1:2) = 0, i.e., liquid-water enthalpy
            HREF_LOC(2) = 1007.0_wp*TREF_LOC
            HREF_LOC(3) = 1870.0_wp*TREF_LOC - 2501600_wp    ! latent heat of vaporization at 273.15 K is 2501.6 kJ/kg

            ! Entropy of Formation in J /kg /K
            SREF_LOC(1) = 0.0_wp
            SREF_LOC(2) = 0.0_wp
            SREF_LOC(3) = -2501600_wp/TREF_LOC

            ! Heat capacity polynomial in J /kg /K; using only the low temperature range
            THERMO_CP(1, :, 1) = 1870.0_wp     ! water vapor
            THERMO_CP(1, :, 2) = 1007.0_wp     ! dry air
            THERMO_CP(1, :, 3) = 4217.6_wp     ! liquid water

            ! -------------------------------------------------------------------
            ! Load thermodynamic data from chemkin file
        case (MIXT_TYPE_CHEMKIN)
            ! call THERMO_READ_CHEMKIN(chemkin_file)
            THERMO_CP = THERMO_CP*RGAS
            NCP = 5             ! Fifth-order polynomial

        end select

        ! 6th and 7th coefficients of heat capacity polynomial are calculated from reference enthalpy
        do is = 1, NSP
            THERMO_CP(6, :, is) = HREF_LOC(is) - THERMO_CP(1, :, is)*TREF_LOC - THERMO_CP(2, :, is)*TREF_LOC*TREF_LOC*0.5_wp
            THERMO_CP(7, :, is) = SREF_LOC(is) - THERMO_CP(2, :, is)*TREF_LOC
        end do

        ! Change heat capacities from molar to mass specific, i.e., J /K /kg
        if (molar_data) then
            do is = 1, NSP
                THERMO_CP(:, :, is) = THERMO_CP(:, :, is)/WGHT(is)
            end do
        end if

        ! ###################################################################
        ! Phase change. Polynomial fit to saturation vapor pressure
        !
        ! p_sat(T) = \sum_1^9 a_i T^{i-1}
        !
        ! ###################################################################
        THERMO_PSAT(:) = 0.0_wp ! Initialize to zero
        NPSAT = 0

        select case (imixture)
            ! Flatau et al., J. Applied Meteorol., 1507-1513, 1992
        case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER)

            NPSAT = 9
            WRK1D_LOC(1) = 0.611213476e+3_wp    ! Pa
            WRK1D_LOC(2) = 0.444007856e+2_wp
            WRK1D_LOC(3) = 0.143064234e+1_wp
            WRK1D_LOC(4) = 0.264461437e-1_wp
            WRK1D_LOC(5) = 0.305930558e-3_wp
            WRK1D_LOC(6) = 0.196237241e-5_wp
            WRK1D_LOC(7) = 0.892344772e-8_wp
            WRK1D_LOC(8) = -0.373208410e-10_wp
            WRK1D_LOC(9) = 0.209339997e-13_wp

            ! going from powers of (T-T_ref) to T, with T_ref = 273.15 K
            TREF_LOC = 273.15_wp
            do ipsat = 1, NPSAT
                THERMO_PSAT(ipsat) = 0.0_wp
                do i = ipsat, NPSAT
                    tmp1 = 1.0_wp
                    do j = i - 1, i - ipsat + 1, -1
                        tmp1 = tmp1*real(j, wp)
                    end do
                    THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat) + &
                                         WRK1D_LOC(i)*TREF_LOC**(i - 1)*tmp1*(-1)**(i - ipsat)
                end do
                tmp2 = 1.0_wp
                do j = ipsat - 1, 1, -1
                    tmp2 = tmp2*real(j, wp)
                end do
                THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/tmp2/TREF_LOC**(ipsat - 1)
            end do

        case default

        end select

        return
    end subroutine Thermo_Mixture_Initialize

    ! ###################################################################
    ! ###################################################################
    ! To be written in terms of elemental functions
    subroutine Thermo_Psat_Polynomial(ijmax, T, p)
        integer(wi) ijmax
        real(wp) T(ijmax)
        real(wp) p(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) i, ipsat

        ! ###################################################################
        if (NPSAT > 0) then
            do i = 1, ijmax
                p(i) = THERMO_PSAT(NPSAT)
                do ipsat = NPSAT - 1, 1, -1
                    p(i) = p(i)*T(i) + THERMO_PSAT(ipsat)
                end do
            end do
        else
            p(:) = 0.0_wp
        end if

        return
    end subroutine Thermo_Psat_Polynomial

    ! ###################################################################
    ! ###################################################################
    subroutine Thermo_dPsat_Polynomial(ijmax, T, dp)
        integer(wi) ijmax
        real(wp) T(ijmax)
        real(wp) dp(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) i, ipsat

        ! ###################################################################
        if (NPSAT > 0) then
            do i = 1, ijmax
                dp(i) = 0.0_wp
                do ipsat = NPSAT - 1, 1, -1
                    dp(i) = dp(i)*T(i) + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                end do
            end do
        else
            dp(:) = 0.0_wp
        end if

        return
    end subroutine Thermo_dPsat_Polynomial

end module Thermo_Base
