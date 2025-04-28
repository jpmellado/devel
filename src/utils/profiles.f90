#include "dns_error.h"

! Definining functions f=f(x) to be used in bcs, ics, and reference background profiles
module Profiles
    use TLab_Constants, only: wp, i4_, pi_wp, efile, wfile, MAX_PARS
    use TLab_WorkFlow,     only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: Profiles_ReadBlock, Profiles_Calculate

    type, public :: profiles_dt
        sequence
        integer type
        integer :: padding = 0_i4_
        logical :: relative = .true.                ! use reference spatial position relative to the extent of the domain
        real(wp) :: mean = 0.0_wp                   ! mean value of f
        real(wp) :: delta = 1.0_wp                  ! increment of f
        real(wp) :: ymean = 0.0_wp                  ! reference spatial position at which f changes      
        real(wp) :: ymean_rel = 0.5_wp              ! same but relative to the extent of the domain
        real(wp) :: thick = 1.0_wp                  ! spatial interval over which f changes
        real(wp) :: lslope = 0.0_wp                 ! slope of f below the ymean
        real(wp) :: uslope = 0.0_wp                 ! slope of f above ymean
        real(wp) :: diam = 0.0_wp                   ! diameter
        real(wp) :: parameters(MAX_PARS) = 0.0_wp   ! additional parameters
    end type profiles_dt

    integer, parameter, public :: PROFILE_NONE               = 0
    integer, parameter, public :: PROFILE_LINEAR             = 1
    integer, parameter, public :: PROFILE_TANH               = 2
    integer, parameter, public :: PROFILE_ERF                = 3
    integer, parameter, public :: PROFILE_BICKLEY            = 4
    integer, parameter, public :: PROFILE_GAUSSIAN           = 5
    integer, parameter, public :: PROFILE_LINEAR_ERF         = 6
    integer, parameter, public :: PROFILE_EKMAN_U            = 7
    integer, parameter, public :: PROFILE_EKMAN_V            = 8
    integer, parameter, public :: PROFILE_EKMAN_U_P          = 9
    integer, parameter, public :: PROFILE_PARABOLIC          = 10
    integer, parameter, public :: PROFILE_LINEAR_CROP        = 11
    integer, parameter, public :: PROFILE_MIXEDLAYER         = 12
    integer, parameter, public :: PROFILE_ERF_ANTISYM        = 13
    integer, parameter, public :: PROFILE_ERF_SURFACE        = 14
    integer, parameter, public :: PROFILE_LINEAR_ERF_SURFACE = 15
    integer, parameter, public :: PROFILE_PARABOLIC_SURFACE  = 16
    integer, parameter, public :: PROFILE_GAUSSIAN_SURFACE   = 17
    integer, parameter, public :: PROFILE_GAUSSIAN_ANTISYM   = 18
    integer, parameter, public :: PROFILE_GAUSSIAN_SYM       = 19
    integer, parameter, public :: PROFILE_TANH_ANTISYM       = 20
    integer, parameter, public :: PROFILE_TANH_SYM           = 21
    integer, parameter, public :: PROFILE_TANH_COS           = 22
    integer, parameter, public :: PROFILE_GAUSSIAN_TANH_SYM  = 23

contains
!########################################################################
!########################################################################
    subroutine Profiles_ReadBlock(bakfile, inifile, block, tag, var, default)
        character(len=*),  intent(in)           :: bakfile, inifile, block, tag
        type(profiles_dt), intent(out)          :: var
        character(len=*),  intent(in), optional :: default

        character(len=512) sRes
        real(wp) derivative

        ! -------------------------------------------------------------------
        call TLab_Write_ASCII(bakfile, '#Profile'//trim(adjustl(tag))//'=<None/Tanh/Erf/Ekman/Parabolic/...>')
        call TLab_Write_ASCII(bakfile, '#'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#YMean'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#YMeanRelative'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#Diam'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#Thick'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#Delta'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#LowerSlope'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#UpperSlope'//trim(adjustl(tag))//'=<value>')

        ! -------------------------------------------------------------------
        if (present(default)) then
            sRes = trim(adjustl(default))
        else
            call ScanFile_Char(bakfile, inifile, block, 'Profile'//trim(adjustl(tag)), 'none', sRes)
        end if
        if (trim(adjustl(sRes)) == 'none') then;                       var%type = PROFILE_NONE
        else if (trim(adjustl(sRes)) == 'tanh') then;                  var%type = PROFILE_TANH
        else if (trim(adjustl(sRes)) == 'tanhsymmetric') then;         var%type = PROFILE_TANH_SYM
        else if (trim(adjustl(sRes)) == 'tanhantisymmetric') then;     var%type = PROFILE_TANH_ANTISYM
        else if (trim(adjustl(sRes)) == 'linear') then;                var%type = PROFILE_LINEAR
        else if (trim(adjustl(sRes)) == 'linearcrop') then;            var%type = PROFILE_LINEAR_CROP
        else if (trim(adjustl(sRes)) == 'erf') then;                   var%type = PROFILE_ERF
        else if (trim(adjustl(sRes)) == 'erfsurface') then;            var%type = PROFILE_ERF_SURFACE
        else if (trim(adjustl(sRes)) == 'erfantisym') then;            var%type = PROFILE_ERF_ANTISYM
        else if (trim(adjustl(sRes)) == 'bickley') then;               var%type = PROFILE_BICKLEY
        else if (trim(adjustl(sRes)) == 'gaussian') then;              var%type = PROFILE_GAUSSIAN
        else if (trim(adjustl(sRes)) == 'gaussiansurface') then;       var%type = PROFILE_GAUSSIAN_SURFACE
        else if (trim(adjustl(sRes)) == 'gaussianvaricose') then;      var%type = PROFILE_GAUSSIAN_ANTISYM
        else if (trim(adjustl(sRes)) == 'gaussiansinuous') then;       var%type = PROFILE_GAUSSIAN_SYM
        else if (trim(adjustl(sRes)) == 'ekman') then;                 var%type = PROFILE_EKMAN_U
        else if (trim(adjustl(sRes)) == 'ekmanp') then;                var%type = PROFILE_EKMAN_U_P
        else if (trim(adjustl(sRes)) == 'parabolic') then;             var%type = PROFILE_PARABOLIC
        else if (trim(adjustl(sRes)) == 'parabolicsurface') then;      var%type = PROFILE_PARABOLIC_SURFACE
        else if (trim(adjustl(sRes)) == 'mixedlayer') then;            var%type = PROFILE_MIXEDLAYER
        else if (trim(adjustl(sRes)) == 'gaussiantanhsymmetric') then; var%type = PROFILE_GAUSSIAN_TANH_SYM
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong '//trim(adjustl(tag))//' profile.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'Mean'//trim(adjustl(tag)), 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') then ! Backwards compatibility
            call ScanFile_Real(bakfile, inifile, block, trim(adjustl(tag)), '0.0', var%mean)
        else
            call ScanFile_Real(bakfile, inifile, block, 'Mean'//trim(adjustl(tag)), '0.0', var%mean)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'YMean'//trim(adjustl(tag)), 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') then
            var%relative = .true.
            call ScanFile_Real(bakfile, inifile, block, 'YMeanRelative'//trim(adjustl(tag)), '0.5', var%ymean_rel)    ! Position in relative coordinates
            ! Backwards compatibility
            call ScanFile_Char(bakfile, inifile, block, 'YCoor'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                call ScanFile_Real(bakfile, inifile, block, 'YCoor'//trim(adjustl(tag)), '0.5', var%ymean_rel)
                call TLab_Write_ASCII(wfile, 'Update tag YCoor to YMeanRelative.')
            end if
        else
            var%relative = .false.
            call ScanFile_Real(bakfile, inifile, block, 'YMean'//trim(adjustl(tag)), '0.0', var%ymean)         ! Position in absolute coordinates
        end if

        call ScanFile_Real(bakfile, inifile, block, 'Thick'//trim(adjustl(tag)), '0.0', var%thick)
        call ScanFile_Real(bakfile, inifile, block, 'Delta'//trim(adjustl(tag)), '0.0', var%delta)
        ! alternative to provide the variable thick in terms of the maximum derivative
        call ScanFile_Char(bakfile, inifile, block, 'Derivative'//trim(adjustl(tag)), 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call ScanFile_Real(bakfile, inifile, block, 'Derivative'//trim(adjustl(tag)), '0.0', derivative)
            call ScanFile_Char(bakfile, inifile, block, 'Thick'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'void') then
                call Profiles_DerToThick(derivative, var)
            end if
            call ScanFile_Char(bakfile, inifile, block, 'Delta'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'void') then
                call Profiles_DerToDelta(derivative, var)
            end if
        end if

        call ScanFile_Real(bakfile, inifile, block, 'LowerSlope'//trim(adjustl(tag)), '0.0', var%lslope)
        call ScanFile_Real(bakfile, inifile, block, 'UpperSlope'//trim(adjustl(tag)), '0.0', var%uslope)
        call ScanFile_Real(bakfile, inifile, block, 'Diam'//trim(adjustl(tag)), '0.0', var%diam)

        call ScanFile_Real(bakfile, inifile, block, 'SurfaceThick'//trim(adjustl(tag)), '1.0', var%parameters(3))
        call ScanFile_Real(bakfile, inifile, block, 'SurfaceDelta'//trim(adjustl(tag)), '0.0', var%parameters(4))
        ! alternative to provide the variable thick in terms of the maximum derivative
        call ScanFile_Char(bakfile, inifile, block, 'SurfaceDerivative'//trim(adjustl(tag)), 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call ScanFile_Real(bakfile, inifile, block, 'SurfaceDerivative'//trim(adjustl(tag)), '0.0', derivative)
            call ScanFile_Char(bakfile, inifile, block, 'SurfaceThick'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'void') then
                call Profiles_DerToThick(derivative, var)
            end if
            call ScanFile_Char(bakfile, inifile, block, 'SurfaceDelta'//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'void') then
                call Profiles_DerToDelta(derivative, var)
            end if
        end if

        return
    end subroutine Profiles_ReadBlock

!########################################################################
!########################################################################
    function Profiles_Calculate(var, y) result(f)
        type(profiles_dt), intent(in) :: var
        real(wp),          intent(in) :: y
        real(wp)                      :: f

        ! -------------------------------------------------------------------
        real(wp) yrel, xi, amplify, zamp, cnought

        ! ###################################################################
        yrel = y - var%ymean    ! position relative to reference height
        amplify = 0.0_wp        ! default

        ! -------------------------------------------------------------------
        ! base state varying between two constant levels
        ! -------------------------------------------------------------------
        if (var%thick == 0.0_wp) then
            if (var%type > 0) amplify = 0.5_wp*sign(1.0_wp, yrel)

        else
            xi = yrel/var%thick

            select case (var%type)

            case (PROFILE_LINEAR)
                amplify = -xi

            case (PROFILE_TANH)
                amplify = 0.5_wp*tanh(-0.5_wp*xi)

            case (PROFILE_TANH_SYM)
                amplify = 0.5_wp*(tanh(-0.5_wp*(xi - 0.5_wp*var%diam/var%thick)) &
                                  + tanh(0.5_wp*(xi + 0.5_wp*var%diam/var%thick)) - 1.0_wp)

            case (PROFILE_TANH_ANTISYM)
                amplify = 0.25_wp*(tanh(-0.5_wp*(xi - 0.5_wp*var%diam/var%thick)) &
                                   - tanh(0.5_wp*(xi + 0.5_wp*var%diam/var%thick)))

            case (PROFILE_ERF, PROFILE_ERF_ANTISYM, PROFILE_ERF_SURFACE)
                amplify = 0.5_wp*erf(-0.5_wp*xi)

            case (PROFILE_PARABOLIC, PROFILE_PARABOLIC_SURFACE)
                amplify = (1.0_wp + 0.5_wp*xi)*(1.0_wp - 0.5_wp*xi)

            case (PROFILE_BICKLEY)
                amplify = 1.0_wp/(cosh(0.5_wp*xi))**2.0_wp

            case (PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SURFACE, PROFILE_GAUSSIAN_TANH_SYM)
                amplify = exp(-0.5_wp*xi**2.0_wp)

            case (PROFILE_GAUSSIAN_SYM)
                amplify = exp(-0.5_wp*(xi - 0.5_wp*var%diam/var%thick)**2.0_wp) &
                          + exp(-0.5_wp*(xi + 0.5_wp*var%diam/var%thick)**2.0_wp)

            case (PROFILE_GAUSSIAN_ANTISYM)
                amplify = exp(-0.5_wp*(xi - 0.5_wp*var%diam/var%thick)**2.0_wp) &
                          - exp(-0.5_wp*(xi + 0.5_wp*var%diam/var%thick)**2.0_wp)

            case (PROFILE_EKMAN_U)
                amplify = 1.0_wp - exp(-xi)*cos(xi)

            case (PROFILE_EKMAN_U_P)
                amplify = 1.0_wp - exp(-xi)*cos(xi) ! + perturbation:

                cnought = pi_wp*pi_wp/4.0_wp/4.0_wp       ! Maximum initial Perturbation is at y=pi/2*var%thick
                zamp = sqrt(2.0_wp)*xi*exp(-xi*xi/8.0_wp/cnought)/(var%thick*var%thick*4.0_wp*cnought)**1.5_wp
                amplify = amplify + zamp                  ! Add Perturbations

            case (PROFILE_EKMAN_V)
                amplify = -exp(-xi)*sin(xi)

            end select

        end if

        ! var%mean profile plus two linear-varying layers
        f = var%mean + var%delta*amplify &
            + var%lslope*yrel*0.5_wp*(1.0_wp - sign(1.0_wp, yrel)) &
            + var%uslope*yrel*0.5_wp*(1.0_wp + sign(1.0_wp, yrel))

        ! -------------------------------------------------------------------
        ! special profiles
        ! -------------------------------------------------------------------
        select case (var%type)

        case (PROFILE_LINEAR_CROP)
            if (yrel < 0.0_wp) then
                f = min(var%lslope*yrel, var%lslope*var%thick)
            else
                f = max(var%uslope*yrel, var%uslope*var%thick)
            end if

        case (PROFILE_MIXEDLAYER)
            if (yrel < 0.0_wp) then
                f = min(var%lslope*yrel, var%lslope*var%thick)
            else
                f = max(var%uslope*yrel, var%uslope*var%thick)
            end if
            f = f - 0.25_wp*var%uslope*var%thick*(1.0_wp - sign(1.0_wp, y - var%thick))

        case (PROFILE_ERF_SURFACE)
            xi = y/var%parameters(3)
            f = f + var%parameters(4)*0.5_wp*(1.0_wp + erf(-0.5_wp*xi))

        case (PROFILE_GAUSSIAN_TANH_SYM)
            amplify = tanh(-0.5_wp*(yrel - 0.5_wp*var%diam)/var%parameters(3)) &
                      + tanh(0.5_wp*(yrel + 0.5_wp*var%diam)/var%parameters(3)) - 1.0_wp
            f = f*amplify

        end select

        return
    end function Profiles_Calculate

!########################################################################
!########################################################################
    subroutine Profiles_DerToThick(derivative, var)  ! Obtain thick from the value of the maximum derivative
        real(wp),          intent(in)    :: derivative
        type(profiles_dt), intent(inout) :: var

        real(wp) thick_ratio    ! for readibility

        select case (var%type)

        case (PROFILE_TANH, PROFILE_TANH_SYM, PROFILE_TANH_ANTISYM)
            thick_ratio = 4.0_wp
            var%thick = -var%delta/derivative/thick_ratio

        case (PROFILE_ERF, PROFILE_ERF_ANTISYM)
            thick_ratio = 2.0_wp*sqrt(pi_wp)
            var%thick = -var%delta/(derivative - var%uslope)/thick_ratio

        case (PROFILE_ERF_SURFACE)
            thick_ratio = 2.0_wp*sqrt(pi_wp)
            var%parameters(3) = -var%parameters(4)/derivative/thick_ratio

        case default
            call TLab_Write_ASCII(efile, __FILE__//'. Undevelop derivative input for this case.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)

        end select

        return
    end subroutine Profiles_DerToThick

!########################################################################
!########################################################################
    subroutine Profiles_DerToDelta(derivative, var)  ! Obtain thick from the value of the maximum derivative
        real(wp),          intent(in)    :: derivative
        type(profiles_dt), intent(inout) :: var

        real(wp) thick_ratio    ! for readibility

        select case (var%type)

        case (PROFILE_TANH, PROFILE_TANH_SYM, PROFILE_TANH_ANTISYM)
            thick_ratio = 4.0_wp
            var%delta = -var%thick*derivative*thick_ratio

        case (PROFILE_ERF, PROFILE_ERF_ANTISYM)
            thick_ratio = 2.0_wp*sqrt(pi_wp)
            var%delta = -var%thick*(derivative - var%uslope)*thick_ratio

        case (PROFILE_ERF_SURFACE)
            thick_ratio = 2.0_wp*sqrt(pi_wp)
            var%parameters(4) = -var%parameters(3)*derivative*thick_ratio

        case default
            call TLab_Write_ASCII(efile, __FILE__//'. Undevelop derivative input for this case.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)

        end select

        return
    end subroutine Profiles_DerToDelta

end module Profiles
