!# Calculate thermodynamic properties of airwater in the anelastic approximation,
!# when a background state is given in the form of appropriate profiles
!#
!# s1 is specific static energy
!# s2 is total water specific humidity
!# s3 is liquid water specific humidity, if any

module Thermo_Anelastic
    use TLab_Constants, only: wp, wi
    use TLab_Memory, only: inb_scal_array
    use Thermo_Base, only: imixture
    use Thermo_Base, only: MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER
    use Thermo_Base, only: THERMO_PSAT, NPSAT, Thermo_Psat_Polynomial
    use Thermo_Base, only: gamma0
    use Thermo_Base, only: nondimensional
    use Thermo_AirWater, only: Rv, Rd, Rdv, Cd, Cdv, Lv0, Ld, Ldv, Cvl, Cdl, Cl, rd_ov_rv, PREF_1000
    use Thermo_AirWater, only: inb_scal_e, inb_scal_ql, inb_scal_T
    implicit none
    private

    public :: Thermo_Anelastic_Initialize
    ! public :: Thermo_Anelastic_Memory

    public :: Thermo_Anelastic_EquilibriumPH
    public :: Thermo_Anelastic_T
    public :: Thermo_Anelastic_StaticL
    public :: Thermo_Anelastic_Buoyancy
    public :: Thermo_Anelastic_Weight_InPlace
    public :: Thermo_Anelastic_Weight_OutPlace
    public :: Thermo_Anelastic_Weight_Add
    public :: Thermo_Anelastic_Weight_Subtract

    public :: Thermo_Anelastic_Rho
    public :: Thermo_Anelastic_Theta
    public :: Thermo_Anelastic_ThetaV
    public :: Thermo_Anelastic_ThetaE
    public :: Thermo_Anelastic_ThetaL
    public :: Thermo_Anelastic_MSE
    public :: Thermo_Anelastic_LapseEquilibrium
    public :: Thermo_Anelastic_LapseFrozen
    public :: Thermo_Anelastic_Pvapor
    public :: Thermo_Anelastic_Weight_DewPoint
    public :: Thermo_Anelastic_RH
    public :: Thermo_Anelastic_STATIC_CONSTANTCP
    ! public :: Thermo_Anelastic_EquilibriumPH_RE

    real(wp), public :: GRATIO = 1.0_wp     ! (gamma0-1)/gamma0 = R0/Cp0
    real(wp), public :: scaleheightinv      ! Normalized gravity, or inverse of pressure scale height.
    !                                       Equivalent to 1/(Fr*RRATIO) in compressible formulation

    ! background, reference profiles
    real(wp), allocatable, public :: epbackground(:)                    ! Potential energy
    real(wp), allocatable, public :: pbackground(:)                     ! Pressure background profile
    real(wp), allocatable, public :: tbackground(:)                     ! Temperature
    real(wp), allocatable, public :: rbackground(:), ribackground(:)    ! Density and its inverse

    ! pointers
    real(wp), pointer, public :: p_hl(:) => null()
    real(wp), pointer, public :: p_qt(:) => null()
    real(wp), pointer, public :: p_ql(:) => null()
    real(wp), pointer, public :: p_T(:) => null()

    ! Equilibrium calculations
    real(wp), public :: NEWTONRAPHSON_ERROR, dsmooth

contains
    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        ! character(len=512) sRes

        !########################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Thermodynamics'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call ScanFile_Real(bakfile, inifile, block, 'ScaleHeight', '0.0', scaleheightinv)
        if (scaleheightinv > 0.0_wp) scaleheightinv = 1.0_wp/scaleheightinv

        !########################################################################
        ! Initialize

        ! Parameters in the evolution equations
        ! Anelastic formulations use p nondimensionalized by reference thermodynamic pressure p_0
        if (nondimensional) then
            GRATIO = (gamma0 - 1.0_wp)/gamma0       ! R_0/C_{p,0}
        end if

        inb_scal_array = inb_scal_array + 1         ! Space for T as diagnostic
        inb_scal_e = 1              ! scalar index for the energy (liquid water static energy)
        select case (imixture)
        case (MIXT_TYPE_AIR)
            inb_scal_T = 2          ! scalar index for the temperature
        case (MIXT_TYPE_AIRVAPOR)
            inb_scal_T = 3          ! scalar index for the temperature
        case (MIXT_TYPE_AIRWATER)
            inb_scal_ql = 3         ! scalar index for the liquid (liquid water specific humidity)
            inb_scal_T = 4          ! scalar index for the temperature
        end select

        return
    end subroutine Thermo_Anelastic_Initialize

    ! !########################################################################
    ! !########################################################################
    ! subroutine Thermo_Anelastic_Memory()
    !     use TLab_Memory, only: isize_field
    !     use TLab_Arrays, only: s

    !     integer(wi) idummy(2)

    !     idummy = shape(s)
    !     ! Pointers
    !     if (idummy(2) >= 1) p_hl(1:isize_field) => s(1:isize_field, 1)
    !     if (idummy(2) >= 2) p_qt(1:isize_field) => s(1:isize_field, 2)
    !     if (idummy(2) >= 3) p_ql(1:isize_field) => s(1:isize_field, 3)
    !     if (idummy(2) >= 4) p_T(1:isize_field) => s(1:isize_field, 4)

    !     ! Shall we add here the creation of the background profiles?

    !     return
    ! end subroutine Thermo_Anelastic_Memory

    !########################################################################
    !########################################################################
    ! Kernel routines. They need to be fast.

    !########################################################################
    !# Calculating the equilibrium T and q_l for given enthalpy and pressure.
    !# Assumes often that THERMO_AI(6,1,1) = THERMO_AI(6,1,2) = 0
    !#
    !# Routine Thermo_Psat_Polynomial is duplicated here to avoid array calls
    !#
    !# Smoothing according to Eq. 25 in Mellado et al., TCFD, 2010
    !#
    !# s1 is total water specific humidity
    !# s2 is liquid water specific humidity
    !#
    !########################################################################
    subroutine Thermo_Anelastic_EquilibriumPH(nx, ny, nz, s, h)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(inout) :: s(nx*ny*nz, 3)   ! qt, ql, T
        real(wp), intent(in) :: h(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) inr, nrmax
        real(wp) ALPHA_1, ALPHA_2, BETA_1, BETA_2, alpha, beta, C_LN2_L
        real(wp) qsat, H_LOC, B_LOC(10), FUN, DER
        real(wp) dsmooth_loc, dqldqt, dqsdt!, qvequ
        real(wp) dummy

        integer(wi) ij, i, k, is, ipsat
        real(wp) psat

        real(wp) E_LOC, P_LOC, T_LOC

        ! ###################################################################
        NEWTONRAPHSON_ERROR = 0.0_wp
        ! maximum number of iterations in Newton-Raphson
        nrmax = 5

        ALPHA_1 = rd_ov_rv*Lv0
        ALPHA_2 = Lv0*(1.0_wp - rd_ov_rv)
        BETA_1 = rd_ov_rv*Cvl + Cd
        BETA_2 = Cdl - rd_ov_rv*Cvl

        C_LN2_L = log(2.0_wp)

        ! ###################################################################
        ij = 0
        do k = 1, nz
            P_LOC = pbackground(k)
            E_LOC = epbackground(k)

            do i = 1, nx*ny
                ij = ij + 1

                H_LOC = h(ij) - E_LOC ! enthalpy

                ! -------------------------------------------------------------------
                ! reference case assuming ql = 0
                ! -------------------------------------------------------------------
                s(ij, 2) = 0.0_wp
                T_LOC = H_LOC/(Cd + s(ij, 1)*Cdv)

                ! calculate saturation specific humidity q_sat(T,p)
                psat = THERMO_PSAT(NPSAT)
                do ipsat = NPSAT - 1, 1, -1
                    psat = psat*T_LOC + THERMO_PSAT(ipsat)
                end do
                dummy = rd_ov_rv/(P_LOC/psat - 1.0_wp)
                qsat = dummy/(1.0_wp + dummy)

                ! -------------------------------------------------------------------
                ! calculate smoothed piecewise linear contribution, if needed
                ! -------------------------------------------------------------------
                if (dsmooth > 0.0_wp) then
                    ! calculate dqsdt from dpsatdt (qs here mean qvequ)
                    dqsdt = 0.0_wp
                    do ipsat = NPSAT - 1, 1, -1
                        dqsdt = dqsdt*T_LOC + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                    end do
                    dqsdt = qsat/psat/(1.0_wp - psat/P_LOC)*dqsdt
                    ! calculate dqldqt
                    dqsdt = dqsdt/(Cd + qsat*Cdv)
                    dqldqt = (1.0_wp/(1.0_wp - qsat) + Cdv*T_LOC*dqsdt)/ &
                             (1.0_wp + (Lv0 - Cvl*T_LOC)*dqsdt)

                    dsmooth_loc = dsmooth*qsat
                    if (s(ij, 1) - qsat < 0.0_wp) then
                        s(ij, 2) = dqldqt*dsmooth_loc*log(exp((s(ij, 1) - qsat)/dsmooth_loc) + 1.0_wp)
                    else
                        s(ij, 2) = dqldqt*((s(ij, 1) - qsat) &
                                           + dsmooth_loc*(C_LN2_L - log(tanh((s(ij, 1) - qsat)/(2.0_wp*dsmooth_loc)) + 1.0_wp)))
                    end if

                end if

                ! -------------------------------------------------------------------
                ! if q_s < q_t, then we have to recalculate T
                ! -------------------------------------------------------------------
                if (qsat < s(ij, 1)) then
                    ! preparing Newton-Raphson
                    alpha = (ALPHA_1 + s(ij, 1)*ALPHA_2 + H_LOC)/P_LOC
                    beta = (BETA_1 + s(ij, 1)*BETA_2)/P_LOC
                    B_LOC(1) = H_LOC + s(ij, 1)*Lv0 - THERMO_PSAT(1)*alpha
                    do is = 2, 9
                        B_LOC(is) = THERMO_PSAT(is - 1)*beta - THERMO_PSAT(is)*alpha
                    end do
                    B_LOC(2) = B_LOC(2) - Cd - s(ij, 1)*Cdl
                    B_LOC(10) = THERMO_PSAT(9)*beta

                    ! executing Newton-Raphson
                    do inr = 1, nrmax
                        FUN = B_LOC(10)
                        DER = 0.0_wp
                        do is = 9, 1, -1
                            FUN = FUN*T_LOC + B_LOC(is)
                            DER = DER*T_LOC + B_LOC(is + 1)*real(is, wp)
                        end do
                        T_LOC = T_LOC - FUN/DER
                    end do
                    NEWTONRAPHSON_ERROR = max(NEWTONRAPHSON_ERROR, abs(FUN/DER)/T_LOC)

                    ! calculate equilibrium vapor specific humidity
                    psat = 0.0_wp
                    do ipsat = NPSAT, 1, -1
                        psat = psat*T_LOC + THERMO_PSAT(ipsat)
                    end do
                    !           dummy = rd_ov_rv /( P_LOC/psat -1.0_wp )
                    !           qvequ = dummy *( 1.0_wp -s(ij,1) )

                    if (dsmooth > 0.0_wp) then ! add correction
                        !              s(ij,2) = s(ij,2) +s(ij,1) -qvequ - (s(ij,1) -qsat) *dqldqt
                        s(ij, 2) = s(ij, 2) + s(ij, 1) - rd_ov_rv/(P_LOC/psat - 1.0_wp)*(1.0_wp - s(ij, 1)) - (s(ij, 1) - qsat)*dqldqt
                    else                           ! or calculate new
                        !              s(ij,2) =          s(ij,1) -qvequ
                        s(ij, 2) = s(ij, 1) - rd_ov_rv/(P_LOC/psat - 1.0_wp)*(1.0_wp - s(ij, 1))
                    end if

                end if

            end do
        end do

        return
    end subroutine Thermo_Anelastic_EquilibriumPH

    !########################################################################
    !########################################################################
#define E_LOC epbackground(k)
#define P_LOC pbackground(k)
#define R_LOC rbackground(k)

    subroutine Thermo_Anelastic_T(nx, ny, nz, s, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz, 3)
        real(wp), intent(out) :: T(nx*ny, nz)

        integer(wi) k

! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_AIR)
            do k = 1, nz
                T(:, k) = s(:, k, 1) - E_LOC
            end do

        case (MIXT_TYPE_AIRVAPOR)
            do k = 1, nz
                T(:, k) = (s(:, k, 1) - E_LOC)/(Cd + s(:, k, 2)*Cdv)
            end do

        case (MIXT_TYPE_AIRWATER)
            do k = 1, nz
                T(:, k) = (s(:, k, 1) - E_LOC + s(:, k, 3)*Lv0)/(Cd + s(:, k, 2)*Cdv + s(:, k, 3)*Cvl)
            end do

        end select

        return
    end subroutine Thermo_Anelastic_T

!########################################################################
! Calculating h_l - h; very similar to the temperature routine
!########################################################################
    subroutine Thermo_Anelastic_StaticL(nx, ny, nz, s, result)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz, *)
        real(wp), intent(out) :: result(nx*ny, nz)

        integer(wi) k

! ###################################################################
#define T_LOC(ij,k) result(ij,k)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            do k = 1, nz
                T_LOC(:, k) = s(:, k, 1) - E_LOC
                result(:, k) = Cl*T_LOC(:, k) + E_LOC - Lv0 - s(:, k, 1)
            end do

        case (MIXT_TYPE_AIRVAPOR)
            do k = 1, nz
                T_LOC(:, k) = (s(:, k, 1) - E_LOC)/(Cd + s(:, k, 2)*Cdv)
                result(:, k) = Cl*T_LOC(:, k) + E_LOC - Lv0 - s(:, k, 1)
            end do

        case (MIXT_TYPE_AIRWATER)
            do k = 1, nz
                T_LOC(:, k) = (s(:, k, 1) - E_LOC + s(:, k, 3)*Lv0)/(Cd + s(:, k, 2)*Cdv + s(:, k, 3)*Cvl)
                result(:, k) = Cl*T_LOC(:, k) + E_LOC - Lv0 - s(:, k, 1)
            end do

        end select

#undef T_LOC

        return
    end subroutine Thermo_Anelastic_StaticL

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_Buoyancy(nx, ny, nz, s, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz, *)
        real(wp), intent(out) :: b(nx*ny, nz)

        integer(wi) k
        real(wp) R_LOC_INV

        ! ###################################################################
#define T_LOC(ij,k) b(ij,k)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            do k = 1, nz
                R_LOC_INV = 1.0_wp/R_LOC
                T_LOC(:, k) = s(:, k, 1) - E_LOC
                b(:, k) = R_LOC_INV*(R_LOC - P_LOC/T_LOC(:, k))
            end do

        case (MIXT_TYPE_AIRVAPOR)
            do k = 1, nz
                R_LOC_INV = 1.0_wp/R_LOC
                T_LOC(:, k) = (s(:, k, 1) - E_LOC)/(Cd + s(:, k, 2)*Cdv)
                T_LOC(:, k) = (Rd + s(:, k, 2)*Rdv)*T_LOC(:, k)                     ! Multiply by gas constant
                b(:, k) = R_LOC_INV*(R_LOC - P_LOC/T_LOC(:, k))
            end do

        case (MIXT_TYPE_AIRWATER)
            do k = 1, nz
                R_LOC_INV = 1.0_wp/R_LOC
                T_LOC(:, k) = (s(:, k, 1) - E_LOC + s(:, k, 3)*Lv0)/(Cd + s(:, k, 2)*Cdv + s(:, k, 3)*Cvl)
                T_LOC(:, k) = (Rd + s(:, k, 2)*Rdv - s(:, k, 3)*Rv)*T_LOC(:, k)     ! Multiply by gas constant
                b(:, k) = R_LOC_INV*(R_LOC - P_LOC/T_LOC(:, k))
            end do

        end select

#undef T_LOC

        return
    end subroutine Thermo_Anelastic_Buoyancy

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_Weight_InPlace(nx, ny, nz, weight, a)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(nz)
        real(wp), intent(inout) :: a(nx*ny, nz)

        integer(wi) k

        ! ###################################################################
        do k = 1, nz
            a(:, k) = a(:, k)*weight(k)
        end do

        return
    end subroutine Thermo_Anelastic_Weight_InPlace

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_Weight_OutPlace(nx, ny, nz, weight, a, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(nz)
        real(wp), intent(in) :: a(nx*ny, nz)
        real(wp), intent(out) :: b(nx*ny, nz)

        integer(wi) k

        ! ###################################################################
        do k = 1, nz
            b(:, k) = a(:, k)*weight(k)
        end do

        return
    end subroutine Thermo_Anelastic_Weight_OutPlace

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_Weight_Add(nx, ny, nz, weight, a, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(nz)
        real(wp), intent(in) :: a(nx*ny, nz)
        real(wp), intent(out) :: b(nx*ny, nz)

        integer(wi) k

        ! ###################################################################
        do k = 1, nz
            b(:, k) = b(:, k) + a(:, k)*weight(k)
        end do

        return
    end subroutine Thermo_Anelastic_Weight_Add

!########################################################################
!########################################################################
    subroutine Thermo_Anelastic_Weight_Subtract(nx, ny, nz, weight, a, b)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: weight(nz)
        real(wp), intent(in) :: a(nx*ny, nz)
        real(wp), intent(out) :: b(nx*ny, nz)

        integer(wi) k

        ! ###################################################################
        do k = 1, nz
            b(:, k) = b(:, k) - a(:, k)*weight(k)
        end do

        return
    end subroutine Thermo_Anelastic_Weight_Subtract

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_P(nx, ny, nz, p)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(out) :: p(nx*ny, nz)

        integer(wi) k

        ! ###################################################################
        do k = 1, nz
            p(:, k) = P_LOC
        end do

        return
    end subroutine Thermo_Anelastic_P

#undef P_LOC
#undef E_LOC
#undef R_LOC

    !########################################################################
    !########################################################################
    ! Procedures for post-processing and initialization.
    ! We could write them explicitly as above using loops and only 1 array.
    ! But we favor clarity over speed.
    ! We write them in terms of a few routines despite requiring more memory and time

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_MSE(nx, ny, nz, s, result)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz, *)
        real(wp), intent(out) :: result(nx*ny, nz)

        integer(wi) k

! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_AIR)
            do k = 1, nz
                result(:, k) = s(:, k, 1)
            end do

        case (MIXT_TYPE_AIRVAPOR)
            do k = 1, nz
                result(:, k) = s(:, k, 1) + s(:, k, 2)*Lv0
            end do

        case (MIXT_TYPE_AIRWATER)
            do k = 1, nz
                result(:, k) = s(:, k, 1) + s(:, k, 2)*Lv0
            end do

        end select

        return
    end subroutine Thermo_Anelastic_MSE

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_Rho(nx, ny, nz, s, rho, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: rho(nx*ny*nz), T(nx*ny*nz)

        ! ###################################################################
#define p rho

        call Thermo_Anelastic_P(nx, ny, nz, p)
        call Thermo_Anelastic_T(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            rho = p/(T*Rd)

        case (MIXT_TYPE_AIRVAPOR)
            rho = p/(T*(Rd + s(:, 2)*Rdv))

        case (MIXT_TYPE_AIRWATER)
            rho = p/(T*(Rd + s(:, 2)*Rdv - s(:, 3)*Rv))

        end select

#undef p

        return
    end subroutine Thermo_Anelastic_Rho

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, pi)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(out) :: pi(nx*ny*nz)

        real(wp) kappa

        ! ###################################################################
        kappa = Rd/Cd*GRATIO

#define p pi

        call Thermo_Anelastic_P(nx, ny, nz, p)

        pi = (PREF_1000/p)**kappa

#undef p

        return
    end subroutine Thermo_Anelastic_ONE_OV_EXNER

    !########################################################################
    ! Dry potential temperature
    !########################################################################
    subroutine Thermo_Anelastic_Theta(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        ! ###################################################################
#define locPi theta

        call Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, locPi)
        call Thermo_Anelastic_T(nx, ny, nz, s, T)

        theta = T*locPi

#undef locPi

        return
    end subroutine Thermo_Anelastic_Theta

    !########################################################################
    ! Virtual Potential temperature
    !########################################################################
    subroutine Thermo_Anelastic_ThetaV(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        real(wp) Rdv_ov_Rd, Rv_ov_Rd

        ! ###################################################################
        Rdv_ov_Rd = Rdv/Rd
        Rv_ov_Rd = Rv/Rd

        call Thermo_Anelastic_Theta(nx, ny, nz, s, theta, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)

        case (MIXT_TYPE_AIRVAPOR)
            theta = theta*(1.0_wp + s(:, 2)*Rdv_ov_Rd)

        case (MIXT_TYPE_AIRWATER)
            theta = theta*(1.0_wp + s(:, 2)*Rdv_ov_Rd - s(:, 3)*Rv_ov_Rd)

        end select

        return
    end subroutine Thermo_Anelastic_ThetaV

    !########################################################################
    ! Liquid water potential temperature
    ! still missing the correction of order one
    !########################################################################
    subroutine Thermo_Anelastic_ThetaL(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        real(wp) Rdv_ov_Rd, Cdv_ov_Cd

        ! ###################################################################
        Rdv_ov_Rd = Rdv/Rd
        Cdv_ov_Cd = Cdv/Cd

#define locPi theta

        call Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, locPi)
        call Thermo_Anelastic_T(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)

        case (MIXT_TYPE_AIRVAPOR)
            theta = T*locPi**((1.0_wp + s(:, 2)*Rdv_ov_Rd)/(1.0_wp + s(:, 2)*Cdv_ov_Cd))

        case (MIXT_TYPE_AIRWATER)
            theta = T*locPi**((1.0_wp + s(:, 2)*Rdv_ov_Rd)/(1.0_wp + s(:, 2)*Cdv_ov_Cd))
            theta = theta*exp(-(Lv0 - T*Cvl)*s(:, 3)/T/(Cd + s(:, 2)*Cdv))

        end select

#undef locPi

        return
    end subroutine Thermo_Anelastic_ThetaL

    !########################################################################
    ! Equivalent potential temperature
    ! still missing the correction of order one
    !########################################################################
    subroutine Thermo_Anelastic_ThetaE(nx, ny, nz, s, theta, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: theta(nx*ny*nz), T(nx*ny*nz)

        real(wp) Cdl_ov_Cd

        ! ###################################################################
        Cdl_ov_Cd = Cdl/Cd

#define locPi theta

        call Thermo_Anelastic_ONE_OV_EXNER(nx, ny, nz, locPi)
        call Thermo_Anelastic_T(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)

        case (MIXT_TYPE_AIRVAPOR)
            theta = T*locPi**((1.0_wp - s(:, 2))/(1.0_wp + s(:, 2)*Cdl_ov_Cd))
            theta = theta*exp((Lv0 - T*Cvl)*s(:, 2)/T/(Cd + s(:, 2)*Cdl))

        case (MIXT_TYPE_AIRWATER)
            theta = T*locPi**((1.0_wp - s(:, 2))/(1.0_wp + s(:, 2)*Cdl_ov_Cd))
            theta = theta*exp((Lv0 - T*Cvl)*(s(:, 2) - s(:, 3))/T/(Cd + s(:, 2)*Cdl))

        end select

#undef locPi

        return
    end subroutine Thermo_Anelastic_ThetaE

    !########################################################################
    !########################################################################
    ! Frozen lapse rate; in unsaturated conditions, this is the unsaturated lapse rate
    subroutine Thermo_Anelastic_LapseFrozen(nx, ny, nz, s, lapse)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: lapse(nx*ny*nz)

        ! ###################################################################
        select case (imixture)
        case (MIXT_TYPE_AIR)
            lapse(:) = GRATIO*scaleheightinv

        case (MIXT_TYPE_AIRVAPOR)
            lapse(:) = GRATIO*scaleheightinv/(Cd + s(:, 2)*Cdv)

        case (MIXT_TYPE_AIRWATER)
            lapse(:) = GRATIO*scaleheightinv/(Cd + s(:, 2)*Cdv + s(:, 3)*Cvl)

        end select

        return
    end subroutine Thermo_Anelastic_LapseFrozen

    !########################################################################
    !########################################################################
    ! Equilibrium lapse rate
    subroutine Thermo_Anelastic_LapseEquilibrium(nx, ny, nz, s, lapse, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: lapse(nx*ny*nz)
        real(wp), intent(inout) :: T(nx*ny*nz)

! -------------------------------------------------------------------
        real(wp) one_ov_Rd, one_ov_Rv, Rv_ov_Rd

! ###################################################################

        select case (imixture)
        case (MIXT_TYPE_AIR)
            call Thermo_Anelastic_LapseFrozen(nx, ny, nz, s, lapse)

        case (MIXT_TYPE_AIRVAPOR)
            call Thermo_Anelastic_LapseFrozen(nx, ny, nz, s, lapse)

        case (MIXT_TYPE_AIRWATER)

#define psat lapse
#define p_ov_psat lapse
#define qv_ov_qd lapse

            call Thermo_Anelastic_T(nx, ny, nz, s, T)

            call Thermo_Psat_Polynomial(nx*ny*nz, T, psat)

            p_ov_psat = 1.0_wp/psat
            call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, pbackground, p_ov_psat)

            qv_ov_qd = rd_ov_rv/(p_ov_psat - 1.0_wp)

            one_ov_Rd = 1.0_wp/(Rd*GRATIO)
            one_ov_Rv = 1.0_wp/(Rv*GRATIO)
            Rv_ov_Rd = Rv/Rd
            lapse = (1.0_wp + qv_ov_qd*(Lv0 - T*Cvl)*one_ov_Rd/T) &
                    /(Cd + s(:, 2)*Cdl - qv_ov_qd*(1.0 - s(:, 2))*Cvl &
                      + qv_ov_qd*(1.0_wp - s(:, 2))*(1.0_wp + qv_ov_qd*Rv_ov_Rd)*(Lv0 - T*Cvl)**2.0_wp*one_ov_Rv/(T*T)) &
                    *GRATIO*scaleheightinv

#undef psat
#undef p_ov_psat
#undef qv_ov_qd

        end select

        return
    end subroutine Thermo_Anelastic_LapseEquilibrium

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_Pvapor(nx, ny, nz, s, pv)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: pv(nx*ny*nz)

        ! ###################################################################
#define p pv
        call Thermo_Anelastic_P(nx, ny, nz, p)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            pv = 0.0_wp

        case (MIXT_TYPE_AIRVAPOR)
            pv = s(:, 2)*Rv/(Rd + s(:, 2)*Rdv)*p

        case (MIXT_TYPE_AIRWATER)
            pv = (s(:, 2) - s(:, 3))*Rv/(Rd + s(:, 2)*Rdv - s(:, 3)*Rv)*p

        end select

#undef p

        return
    end subroutine Thermo_Anelastic_Pvapor

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_RH(nx, ny, nz, s, rh, T)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: rh(nx*ny*nz)
        real(wp), intent(inout) :: T(nx*ny*nz)

        ! ###################################################################
#define pv T
#define psat rh
        call Thermo_Anelastic_T(nx, ny, nz, s, T)
        call Thermo_Psat_Polynomial(nx*ny*nz, T, psat)

        call Thermo_Anelastic_Pvapor(nx, ny, nz, s, pv)

        rh = pv/psat*100.0_wp

#undef pv
#undef psat

        return
    end subroutine Thermo_Anelastic_RH

    !########################################################################
    !########################################################################
    subroutine Thermo_Anelastic_Weight_DewPoint(nx, ny, nz, s, dpvdy, Td, Lapse)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *), dpvdy(nx*ny*nz)
        real(wp), intent(out) :: lapse(nx*ny*nz), Td(nx*ny*nz)

        ! -------------------------------------------------------------------
        integer(wi) inr, nrmax, ij, ipsat
        real(wp) dpsat, psat

        ! ###################################################################
        if (imixture == MIXT_TYPE_AIR) return

        ! maximum number of iterations in Newton-Raphson
        nrmax = 5

#define pv Lapse

        call Thermo_Anelastic_Pvapor(nx, ny, nz, s, pv)

        call Thermo_Anelastic_T(nx, ny, nz, s, Td)    ! Using actual temperature as initial condition

        ! If liquid, we should use qsat instead of ql (s(:,3)) because the smoothing function imposes
        ! an exponentially small value, but nonzero (?)

        do ij = 1, nx*ny*nz
            do inr = 1, nrmax                                   ! executing Newton-Raphson
                psat = THERMO_PSAT(NPSAT); dpsat = 0.0_wp
                do ipsat = NPSAT - 1, 1, -1
                    psat = psat*Td(ij) + THERMO_PSAT(ipsat)
                    dpsat = dpsat*Td(ij) + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                end do
                Td(ij) = Td(ij) - (psat - pv(ij))/dpsat
            end do

            Lapse(ij) = -dpvdy(ij)/dpsat

        end do

#undef pv

        return
    end subroutine Thermo_Anelastic_Weight_DewPoint

    !########################################################################
    !########################################################################
    ! Just to check what the effect of using a wrong cp would be
    subroutine Thermo_Anelastic_STATIC_CONSTANTCP(nx, ny, nz, s, result)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: result(nx*ny*nz)

        ! ###################################################################
#define T result
        call Thermo_Anelastic_T(nx, ny, nz, s, T)

        select case (imixture)
        case (MIXT_TYPE_AIR)
            result = s(:, 1)

        case (MIXT_TYPE_AIRVAPOR)
            result = s(:, 1) + (Cd - (Cd + s(:, 2)*Cdv))*T

        case (MIXT_TYPE_AIRWATER)
            result = s(:, 1) + (Cd - (Cd + s(:, 2)*Cdv + s(:, 3)*Cvl))*T

        end select

#undef T

        return
    end subroutine Thermo_Anelastic_STATIC_CONSTANTCP

end module Thermo_Anelastic
