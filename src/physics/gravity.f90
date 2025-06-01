#include "tlab_error.h"

! Incompressible formulation uses buoyancy and different forms of buoyancy function
! Anelastic formulation uses (rho-rho_ref) *g and scaleheight defines rho_ref
! Compressible formulation uses simply the gravity force rho *g

module Gravity
    use TLab_Constants, only: wp, wi, small_wp, efile, lfile, wfile, MAX_PROF, MAX_VARS, MAX_PARS
    use TLab_Memory, only: inb_scal, inb_scal_array, inb_flow, inb_flow_array
    use NavierStokes, only: froude
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    type gravity_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) vector(3)
    end type gravity_dt
    type(gravity_dt), public, protected :: gravityProps
    real(wp), allocatable, public :: bbackground(:)

    integer, parameter :: TYPE_GRAV_NONE = 0
    integer, parameter :: TYPE_GRAV_HOMOGENEOUS = 5
    integer, parameter :: TYPE_GRAV_LINEAR = 6
    integer, parameter :: TYPE_GRAV_BILINEAR = 7
    integer, parameter :: TYPE_GRAV_QUADRATIC = 8

    public :: Gravity_Initialize
    public :: Gravity_Source
    public :: Gravity_Hydrostatic_Enthalpy

contains
    !########################################################################
    !########################################################################
    subroutine Gravity_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer(wi) idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Gravity'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/Homogeneous/Linear/Bilinear/Quadratic>')
        call TLab_Write_ASCII(bakfile, '#Vector=<Gx,Gy,Gz>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; gravityProps%type = TYPE_GRAV_NONE
        else if (trim(adjustl(sRes)) == 'homogeneous') then; gravityProps%type = TYPE_GRAV_HOMOGENEOUS
        else if (trim(adjustl(sRes)) == 'linear') then; gravityProps%type = TYPE_GRAV_LINEAR
        else if (trim(adjustl(sRes)) == 'bilinear') then; gravityProps%type = TYPE_GRAV_BILINEAR
        else if (trim(adjustl(sRes)) == 'quadratic') then; gravityProps%type = TYPE_GRAV_QUADRATIC
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TermBodyForce option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        gravityProps%vector = 0.0_wp
        call ScanFile_Char(bakfile, inifile, block, 'Vector', '0.0,0.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, gravityProps%vector)

        gravityProps%active = .false.
        if (abs(gravityProps%vector(1)) > 0.0_wp) then; gravityProps%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Gravity along Ox.'); end if
        if (abs(gravityProps%vector(2)) > 0.0_wp) then; gravityProps%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Gravity along Oy.'); end if
        if (abs(gravityProps%vector(3)) > 0.0_wp) then; gravityProps%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Gravity along Oz.'); end if

        if (froude > 0.0_wp) then   ! adding the froude number into the vector g
            gravityProps%vector(:) = gravityProps%vector(:)/froude
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Froude number must be nonzero.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (gravityProps%type /= TYPE_GRAV_NONE) then
            gravityProps%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '0.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, gravityProps%parameters)
            gravityProps%scalar(1) = idummy                                     ! number of scalars affecting buoyancy function
            gravityProps%scalar(1) = min(inb_scal_array, gravityProps%scalar(1))

        end if

        return
    end subroutine Gravity_Initialize

    !########################################################################
    !########################################################################
    subroutine Gravity_Source(locProps, nx, ny, nz, s, b)
        type(gravity_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal_array)
        real(wp), intent(out) :: b(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) k, is
        real(wp) c0_loc, c1_loc, c2_loc, dummy

#define ref(k) bbackground(k)

        ! #######################################################################
        select case (locProps%type)
        case (TYPE_GRAV_HOMOGENEOUS)
            b(:, :, :) = locProps%parameters(1)

        case (TYPE_GRAV_LINEAR)
            c1_loc = locProps%parameters(1)                     ! proportionality factors
            c2_loc = locProps%parameters(2)
            c0_loc = locProps%parameters(inb_scal_array + 1)    ! independent term

            select case (locProps%scalar(1))
            case (1)
                do k = 1, nz
                    dummy = ref(k) - c0_loc
                    b(:, :, k) = c1_loc*s(:, :, k, 1) - dummy
                end do

            case (2)
                do k = 1, nz
                    dummy = ref(k) - c0_loc
                    b(:, :, k) = c1_loc*s(:, :, k, 1) + c2_loc*s(:, :, k, 2) - dummy
                end do

            case default
                do k = 1, nz
                    b(:, :, k) = c0_loc - ref(k)
                    do is = 1, locProps%scalar(1)
                        if (abs(locProps%parameters(is)) > small_wp) b(:, :, k) = b(:, :, k) + locProps%parameters(is)*s(:, :, k, is)
                    end do
                end do

            end select

        case (TYPE_GRAV_BILINEAR)
            c0_loc = locProps%parameters(1); 
            c1_loc = locProps%parameters(2); 
            c2_loc = locProps%parameters(3)

            do k = 1, nz
                b(:, :, k) = c0_loc*s(:, :, k, 1) + c1_loc*s(:, :, k, 2) + c2_loc*s(:, :, k, 1)*s(:, :, k, 2) - ref(k)
            end do

        case (TYPE_GRAV_QUADRATIC)
            c0_loc = -locProps%parameters(1)/(locProps%parameters(2)/2.0_wp)**2
            c1_loc = locProps%parameters(2)

            do k = 1, nz
                b(:, :, k) = c0_loc*s(:, :, k, 1)*(s(:, :, k, 1) - c1_loc) - ref(k)
            end do

        case default
            b = 0.0_wp

        end select

#undef ref

        return
    end subroutine Gravity_Source

    !########################################################################
    ! Compute hydrostatic equilibrium from profiles s=(h,q_t) where h is the enthalpy
    ! Evaluate the integral \int_zref^z dx/H(x), where H(x) is the scale height in the system
    !########################################################################
    subroutine Gravity_Hydrostatic_Enthalpy(fdmi, z, s, ep, T, p, zref, pref, wrk1d)
        use TLab_Constants, only: BCS_MIN
        use FDM_Integral, only: FDM_Int1_Solve, fdm_integral_dt
        use Thermodynamics, only: imode_thermo, THERMO_TYPE_ANELASTIC, THERMO_TYPE_COMPRESSIBLE
        use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
        use Thermo_Anelastic
        use Microphysics
        use OPR_ODES

        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: z(:)
        real(wp), intent(inout) :: s(:, :)      ! We calculate equilibrium composition
        real(wp), intent(out) :: ep(:), T(:), p(:)
        real(wp), intent(in) :: zref, pref
        real(wp), dimension(size(z), 3), intent(inout) :: wrk1d

        ! -------------------------------------------------------------------
        integer(wi) iter, niter, k, kcenter, nz
        real(wp) dummy

        ! ###################################################################
        nz = size(z)

        ! Get the center
        do k = 1, nz
            if (z(k) <= zref .and. &
                z(k + 1) > zref) then
                kcenter = k
                exit
            end if
        end do

        ! specific potential energy
        if (imode_thermo == THERMO_TYPE_ANELASTIC) then
            ep(:) = (z(:) - zref)*GRATIO*scaleheightinv
            epbackground(:) = ep(:)
        else
            ep(:) = -(z(:) - zref)*gravityProps%vector(3)
        end if

        ! hydrstatic pressure
#define p_aux(i)        wrk1d(i,1)
#define r_aux(i)        wrk1d(i,2)
#define wrk_aux(i)      wrk1d(i,3)

        ! Setting the pressure entry to 1 to get 1/RT
        p_aux(:) = 1.0_wp

        niter = 10

        p(:) = pref                                                                     ! initialize iteration
        s(:, inb_scal + 1:inb_scal_array) = 0.0_wp                                      ! initialize diagnostic
        do iter = 1, niter           ! iterate
            if (imode_thermo == THERMO_TYPE_ANELASTIC) then
                pbackground(:) = p_aux(:)
                call Thermo_Anelastic_Rho(1, 1, nz, s, r_aux(:), wrk_aux(:))    ! Get r_aux=1/RT
                r_aux(:) = -scaleheightinv*r_aux(:)
            else
                ! call THERMO_AIRWATER_PH_RE(nx, s(:, 2), p, s(:, 1), T)
                ! call THERMO_THERMAL_DENSITY(nx, s(:, 2), p_aux(:), T, r_aux(:))     ! Get r_aux=1/RT
                ! r_aux(:) = gravityProps%vector(2)*r_aux(:)
            end if

            p(1) = 0.0_wp
            call FDM_Int1_Solve(1, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, r_aux(:), p, wrk_aux(:))

            ! Calculate pressure and normalize s.t. p=pref at y=zref
            p(:) = exp(p(:))
            if (abs(zref - z(kcenter)) == 0.0_wp) then
                dummy = p(kcenter)
            else
                dummy = p(kcenter) + (p(kcenter + 1) - p(kcenter)) &
                        /(z(kcenter + 1) - z(kcenter))*(zref - z(kcenter))
            end if
            dummy = pref/dummy
            p(:) = dummy*p(:)

            select case (imode_thermo)
            case (THERMO_TYPE_ANELASTIC)
                if (imixture == MIXT_TYPE_AIRWATER .and. evaporationProps%type == TYPE_EVA_EQUILIBRIUM) then
                    pbackground(:) = p(:)
                    call Thermo_Anelastic_EquilibriumPH(1, 1, nz, s(:, 2), s(:, 1))
                end if
                call Thermo_Anelastic_T(1, 1, nz, s, T)

            case (THERMO_TYPE_COMPRESSIBLE)

            end select

        end do

#undef p_aux
#undef r_aux

        return
    end subroutine Gravity_Hydrostatic_Enthalpy

end module Gravity
