#include "tlab_error.h"

! Evolution (prognostic) equations:
! q(:,1:inb_flow) for flow variables.
! s(:,1:inb_scal) for scal variables.

module NavierStokes     ! Shall we call it EvolutionEquations?
    use TLab_Constants, only: wp, wi, lfile, efile, wfile, MAX_VARS
    implicit none
    private

    public :: NavierStokes_Initialize_Parameters

    integer, public, protected :: nse_eqns                                  ! formulation of evolution equations
    integer, parameter, public :: DNS_EQNS_BOUSSINESQ = 1
    integer, parameter, public :: DNS_EQNS_ANELASTIC = 2
    integer, parameter, public :: DNS_EQNS_COMPRESSIBLE = 3

    integer, public, protected :: nse_advection, nse_viscous, nse_diffusion ! formulation of Burgers operator
    integer, parameter, public :: EQNS_NONE = 0
    integer, parameter, public :: EQNS_DIVERGENCE = 1
    integer, parameter, public :: EQNS_SKEWSYMMETRIC = 2
    integer, parameter, public :: EQNS_CONVECTIVE = 3
    integer, parameter, public :: EQNS_EXPLICIT = 4

    ! Nondimensional numbers
    real(wp), public :: visc, schmidt(MAX_VARS)                     ! molecular transport
    real(wp), public, protected :: froude                           ! gravity force
    real(wp), public, protected :: rossby                           ! Coriolis force

    real(wp), public, protected :: prandtl                          ! molecular transport
    real(wp), public, protected :: mach                             ! Mach number

    real(wp), public, protected :: damkohler(MAX_VARS)              ! reaction

    real(wp), public, protected :: stokes                           ! particle inertial effects
    real(wp), public, protected :: settling                         ! sedimentation effects

contains
    subroutine NavierStokes_Initialize_Parameters(inifile)
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        use TLab_Memory, only: inb_flow, inb_flow_array, inb_scal, inb_scal_array
        use TLab_Memory, only: inb_wrk1d, inb_wrk2d
        ! use Thermodynamics, only: mach

        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        character(len=64) lstr
        integer(wi) is, idummy
        real(wp) dummy, reynolds

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'EvolutionEquations'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Equations=<total/internal/incompressible/anelastic>')
        call TLab_Write_ASCII(bakfile, '#TermAdvection=<divergence/skewsymmetric>')
        call TLab_Write_ASCII(bakfile, '#TermViscous=<divergence/explicit>')
        call TLab_Write_ASCII(bakfile, '#TermDiffusion=<divergence/explicit>')

        call ScanFile_Char(bakfile, inifile, block, 'Equations', 'internal', sRes)
        if (trim(adjustl(sRes)) == 'boussinesq') then; nse_eqns = DNS_EQNS_BOUSSINESQ
        else if (trim(adjustl(sRes)) == 'anelastic') then; nse_eqns = DNS_EQNS_ANELASTIC
        else if (trim(adjustl(sRes)) == 'compressible') then; nse_eqns = DNS_EQNS_COMPRESSIBLE
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong entry Main.Equations option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TermAdvection', 'convective', sRes)
        if (trim(adjustl(sRes)) == 'none') then; nse_advection = EQNS_NONE
        else if (trim(adjustl(sRes)) == 'divergence') then; nse_advection = EQNS_DIVERGENCE
        else if (trim(adjustl(sRes)) == 'skewsymmetric') then; nse_advection = EQNS_SKEWSYMMETRIC
        else if (trim(adjustl(sRes)) == 'convective') then; nse_advection = EQNS_CONVECTIVE
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TermAdvection option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TermViscous', 'explicit', sRes)
        if (trim(adjustl(sRes)) == 'none') then; nse_viscous = EQNS_NONE
        else if (trim(adjustl(sRes)) == 'divergence') then; nse_viscous = EQNS_DIVERGENCE
        else if (trim(adjustl(sRes)) == 'explicit') then; nse_viscous = EQNS_EXPLICIT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TermViscous option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TermDiffusion', 'explicit', sRes)
        if (trim(adjustl(sRes)) == 'none') then; nse_diffusion = EQNS_NONE
        else if (trim(adjustl(sRes)) == 'divergence') then; nse_diffusion = EQNS_DIVERGENCE
        else if (trim(adjustl(sRes)) == 'explicit') then; nse_diffusion = EQNS_EXPLICIT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TermDiffusion option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! -------------------------------------------------------------------
        call TLab_Write_ASCII(bakfile, '#Reynolds=<value>')
        call TLab_Write_ASCII(bakfile, '#Schmidt=<value>')
        call TLab_Write_ASCII(bakfile, '#Froude=<value>')
        call TLab_Write_ASCII(bakfile, '#Rossby=<value>')
        call TLab_Write_ASCII(bakfile, '#Damkohler=<value>')
        call TLab_Write_ASCII(bakfile, '#Stokes=<value>')
        call TLab_Write_ASCII(bakfile, '#Settling=<value>')
        call TLab_Write_ASCII(bakfile, '#Mach=<value>')
        call TLab_Write_ASCII(bakfile, '#Prandtl=<value>')

        ! Molecular transport
        call ScanFile_Real(bakfile, inifile, block,  'Reynolds', '-1.0', reynolds)
        if (reynolds <= 0.0) then
            call ScanFile_Real(bakfile, inifile, block,  'Viscosity', '-1.0', dummy)
            if (dummy <= 0.0) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Molecular transport coefficients need to be positive.')
                call TLab_Stop(DNS_ERROR_OPTION)
            else
                reynolds = 1.0_wp/dummy
            end if
        end if

        if (nse_viscous == EQNS_NONE) then
            visc = 0.0_wp
        else
            visc = 1.0_wp/reynolds
        end if

        call ScanFile_Char(bakfile, inifile, block, 'Schmidt', '1.0', sRes)
        schmidt(:) = 0.0_wp; inb_scal = MAX_VARS
        call LIST_REAL(sRes, inb_scal, schmidt)

        ! Gravity
        call ScanFile_Real(bakfile, inifile, block, 'Froude', '-1.0', froude)
        if (froude <= 0.0) then
            call ScanFile_Real(bakfile, inifile, block, 'Gravity', '1.0', dummy)   ! default value
            froude = 1.0_wp/dummy
        end if

        ! Coriolis
        call ScanFile_Real(bakfile, inifile, block, 'Rossby', '-1.0', rossby)
        if (rossby <= 0.0) then
            call ScanFile_Real(bakfile, inifile, block, 'Coriolis', '1.0', dummy)   ! default value
            rossby = 1.0_wp/dummy
        end if

        ! Additional source terms; this can be used to control input in chemistry, radiation, microphysics.... Still needed?
        lstr = '0.0'
        do is = 2, inb_scal
            lstr = trim(adjustl(lstr))//',0.0'
        end do
        call ScanFile_Char(bakfile, inifile, block, 'Damkohler', lstr, sRes)
        damkohler(:) = 0.0_wp; idummy = MAX_VARS
        call LIST_REAL(sRes, idummy, damkohler)
        if (inb_scal /= idummy) then ! Consistency check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Schmidt and Damkholer sizes do not match.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! Compressible flows
        call ScanFile_Real(bakfile, inifile, block, 'Prandtl', '1.0', prandtl)   ! molecular transport, but only appearing in compressible formulation
        call ScanFile_Real(bakfile, inifile, block, 'Mach', '1.0', mach)

        ! Particle-laden flows
        call ScanFile_Real(bakfile, inifile, block, 'Stokes', '0.0', stokes)
        call ScanFile_Real(bakfile, inifile, block, 'Settling', '0.0', settling)

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        ! prognostic and diagnostic variables
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            inb_flow = 3                            ! space for u, v, w
            inb_flow_array = inb_flow

        case (DNS_EQNS_COMPRESSIBLE)
            inb_flow = 5                            ! space for u, v, w, e, rho
            inb_flow_array = inb_flow + 2           ! space for p, T

        end select

        inb_scal_array = inb_scal ! Default is that q/s arrays contain only the prognostic variables;
        !                           can be changed in Thermo_Initialize(ifile)

        ! scratch arrays
        inb_wrk1d = 18

        inb_wrk2d = 3

        return
    end subroutine NavierStokes_Initialize_Parameters

end module NavierStokes
