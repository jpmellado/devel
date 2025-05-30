#include "tlab_error.h"

module Radiation
    use TLab_Constants, only: wp, wi, pi_wp, efile, MAX_PARS, MAX_VARS
    use TLab_Constants, only: BCS_MAX, BCS_MIN
    use TLab_Grid, only: z
    use FDM_Integral, only: FDM_Int1_Solve, fdm_integral_dt
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    ! use TLab_Memory, only: inb_scal_array
    use TLab_Arrays, only: wrk2d, wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_AirWater, only: inb_scal_ql, inb_scal_T
    use OPR_ODES
    use Integration
    implicit none
    private

    type, public :: radterm_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type radterm_dt

    integer, parameter :: ncomps_max = 10                ! maximum number of radiatively active components
    integer, parameter :: nbands_max = 10                ! maximum number of spectral bands
    ! type infrared_dt
    !     sequence
    !     integer type
    !     logical active(MAX_VARS), lpadding(3)               ! evolution equations affected by this term
    !     integer ncomps                                      ! number of radiatively active components
    !     integer nbands                                      ! number of spectral bands
    !     real(wp) :: kappa(ncomps_max, nbands_max)           ! mass absorption coefficients for each radiatively active component
    !     real(wp) beta(3, nbands_max)                        ! polynomial coefficients for band emission functions; assuming second-order polynomial
    !     real(wp) bcs_t(3, nbands_max)                       ! downward fluxes at the top of the domain
    !     real(wp) bcs_b(3, nbands_max)                       ! upward fluxes at the bottom of the domain
    !     real(wp) :: epsilon_b                               ! surface emissivity at ymin
    ! end type infrared_dt

    type(radterm_dt), public, protected :: infraredProps               ! Radiation parameters
    ! type(radterm_dt), public, protected :: visibleProps                ! Radiation parameters

    public :: Radiation_Initialize
    public :: Radiation_Infrared_Z

    integer, parameter, public :: TYPE_RAD_NONE = 0
    integer, parameter :: TYPE_IR_GRAY_LIQUID = 1
    integer, parameter :: TYPE_IR_GRAY = 2
    integer, parameter :: TYPE_IR_BAND = 3
    integer, parameter :: TYPE_BULK1DLOCAL = 10         ! backwards compatibility, to be removed

    real(wp), parameter :: sigma = 5.67037442e-8_wp     ! Stefan-Boltzmann constant, W /m^2 /K^4
    real(wp) :: mu                                      ! mean direction parameter
    ! to be moved into module
    real(wp) :: epsilon                                 ! surface emissivity at ymin
    integer :: ncomps                                   ! number of radiatively active components
    integer :: nbands                                   ! number of spectral bands
    real(wp) beta(3, nbands_max)                        ! polynomial coefficients for band emission functions; assuming second-order polynomial
    real(wp) kappa(ncomps_max, nbands_max)              ! mass absorption coefficients for each radiatively active component in each band
    !
    real(wp), allocatable, target :: bcs_ht(:)          ! flux boundary condition at the top of the domain
    real(wp), allocatable, target :: bcs_hb(:)          ! flux boundary condition at the bottom of the domain
    ! real(wp), allocatable, target :: t_ht(:)            ! temperature at the top of the domain
    real(wp), allocatable, target :: tmp_rad(:, :)      ! 3D temporary arrays for radiation routine

    real(wp), pointer :: p_tau(:, :) => null()
    real(wp), pointer :: p_source(:) => null(), p_flux_up(:) => null(), p_flux_down(:) => null()

contains
    !########################################################################
    !########################################################################
    subroutine Radiation_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, isize_field
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer(wi) idummy, iband, ic
        integer(wi) :: inb_tmp_rad = 0
        real(wp) :: dummy(MAX_PARS)

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Infrared'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Scalar=<value>')
        call TLab_Write_ASCII(bakfile, '#AbsorptionComponent#=<values>')
        call TLab_Write_ASCII(bakfile, '#BoundaryConditions=<values>')
        call TLab_Write_ASCII(bakfile, '#BetaCoefficient=<values>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; infraredProps%type = TYPE_RAD_NONE
        else if (trim(adjustl(sRes)) == 'grayliquid') then; infraredProps%type = TYPE_IR_GRAY_LIQUID
        else if (trim(adjustl(sRes)) == 'gray') then; infraredProps%type = TYPE_IR_GRAY
        else if (trim(adjustl(sRes)) == 'band') then; infraredProps%type = TYPE_IR_BAND
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        infraredProps%active = .false.
        if (infraredProps%type /= TYPE_RAD_NONE) then
            call ScanFile_Int(bakfile, inifile, block, 'Scalar', '1', idummy)         ! in which evolution equation radiation acts
            infraredProps%active(idummy) = .true.

            infraredProps%auxiliar(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'BoundaryConditions', '1.0, 1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, infraredProps%auxiliar)
            epsilon = infraredProps%auxiliar(idummy)        ! last value is surface emissivity at ymin
            nbands = idummy - 1

            kappa(:, :) = 0.0_wp
            do ncomps = 1, ncomps_max
                write (sRes, *) ncomps
                call ScanFile_Char(bakfile, inifile, block, 'AbsorptionComponent'//trim(adjustl(sRes)), 'void', sRes)
                if (trim(adjustl(sRes)) /= 'void') then
                    idummy = nbands_max
                    call LIST_REAL(sRes, idummy, dummy)
                    if (idummy /= nbands) then
                        call TLab_Write_ASCII(efile, __FILE__//'. Error in '//trim(adjustl(block))//'.AbsorptionComponent.')
                        call TLab_Stop(DNS_ERROR_OPTION)
                    end if
                    kappa(ncomps, 1:nbands) = dummy(1:nbands)
                else
                    exit
                end if
            end do
            ncomps = ncomps - 1           ! correct for the increment in the loop

            beta(:, :) = 0.0_wp
            beta(1:3, 1) = [2.6774e-1_wp, -1.3344e-3_wp, 1.8017e-6_wp]  ! default coefficients for band 1 according to Jeevanjee, 2023 for vapor bands
            beta(1:3, 2) = [-2.2993e-2_wp, 8.7439e-5_wp, 1.4744e-7_wp]  ! default coefficients for band 2
            do ic = 1, 3                                                ! so far, only 3 coefficients, 2. order polynomial
                write (sRes, *) ic
                call ScanFile_Char(bakfile, inifile, block, 'BetaCoefficient'//trim(adjustl(sRes)), 'void', sRes)
                if (trim(adjustl(sRes)) /= 'void') then
                    idummy = nbands_max
                    call LIST_REAL(sRes, idummy, dummy)
                    if (idummy /= nbands - 1) then
                        call TLab_Write_ASCII(efile, __FILE__//'. Error in '//trim(adjustl(block))//'.BetaCoefficient.')
                        call TLab_Stop(DNS_ERROR_OPTION)
                    end if
                    beta(ic, 1:nbands) = dummy(1:nbands)
                end if
            end do
            beta(1:3, nbands) = [1.0_wp, 0.0_wp, 0.0_wp]                ! last band from equation sum beta_i = 1
            do iband = 1, nbands - 1
                beta(1:3, nbands) = beta(1:3, nbands) - beta(1:3, iband)
            end do

            ! -------------------------------------------------------------------
            ! infraredProps%scalar = inb_scal_array                       ! By default, radiation is caused by last scalar

            select case (imixture)
            case (MIXT_TYPE_AIRWATER)
                infraredProps%active(1) = .true.           ! energy
                ! infraredProps%active(inb_scal_h) = .true.           ! energy
                infraredProps%active(inb_scal_ql) = .true.           ! liquid
                infraredProps%active(inb_scal_T) = .true.            ! liquid
                ! first radiatively active scalar is liquid
                ! second radiatively active scalar is vapor
                ! third radiatively active scalar is a homogeneous field, e.g., CO2

            case default
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Infrared only derived for airwater mixture.')
                call TLab_Stop(DNS_ERROR_OPTION)

            end select

        end if

        mu = 0.5_wp*(1.0_wp/sqrt(3.0_wp) + 1.0_wp/sqrt(2.0_wp))         ! mean direction, in (1/sqrt{3}, 1/sqrt{2})
        ! mu = 1.0_wp/sqrt(2.0_wp)
        ! mu = 0.5_wp     ! testing

        ! !########################################################################
        ! bakfile = trim(adjustl(inifile))//'.bak'
        ! block = 'Visible'

        !########################################################################
        ! Local allocation
        allocate (bcs_ht(imax*jmax), bcs_hb(imax*jmax))!, t_ht(imax*jmax))
        select case (infraredProps%type)
        case (TYPE_IR_BAND)
            inb_tmp_rad = 3                                         ! Additional memory space
            call TLab_Allocate_Real(__FILE__, tmp_rad, [isize_field, inb_tmp_rad], 'tmp-rad')
            p_source(1:isize_field) => tmp_rad(1:isize_field, 1)
            p_flux_up(1:isize_field) => tmp_rad(1:isize_field, 2)
            p_flux_down(1:isize_field) => tmp_rad(1:isize_field, 3)
        end select

        return
    end subroutine Radiation_Initialize

    !########################################################################
    !########################################################################
    subroutine Radiation_Infrared_Z(localProps, nx, ny, nz, fdmi, s, source, b, flux_down, flux_up, flux)
        use Thermo_Anelastic
        type(radterm_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nx, ny, nz
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: source(nx*ny*nz)           ! also used for absorption coefficient
        real(wp), intent(inout) :: b(nx*ny*nz)              ! emission function, returns flux up
        real(wp), intent(inout) :: flux_down(nx*ny*nz), flux_up(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)   ! net flux upwards

        ! target :: source, b, tmp1, flux

        ! -----------------------------------------------------------------------
        integer iband

        !########################################################################
        p_tau(1:nx*ny, 1:nz) => wrk3d(1:nx*ny*ny)                 ! set pointer to optical depth and transmission functions used in routines below

        ! -----------------------------------------------------------------------
        select case (localProps%type)
        case (TYPE_IR_GRAY_LIQUID)
            source(1:nx*ny*nz) = kappa(1, 1)*s(:, inb_scal_ql)  ! absorption coefficient in array source to save memory
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, source)
            end if

            bcs_ht = localProps%auxiliar(1)                     ! downward flux at domain top
            bcs_hb = localProps%auxiliar(2)                     ! upward flux at domain bottom
            if (present(flux)) then                             ! solve radiative transfer equation along y
                call IR_RTE1_OnlyLiquid(localProps, nx*ny, nz, fdmi, source, flux_down, flux_up)
            else
                call IR_RTE1_OnlyLiquid(localProps, nx*ny, nz, fdmi, source)
            end if

            ! -----------------------------------------------------------------------
        case (TYPE_IR_GRAY)
            b(1:nx*ny*nz) = sigma*s(1:nx*ny*nz, inb_scal_T)**4.0_wp             ! emission function, Stefan-Boltzmann law

            source(1:nx*ny*nz) = kappa(1, 1)*s(:, inb_scal_ql) + &              ! absorption coefficient
                                 kappa(2, 1)*(s(:, 2) - s(:, inb_scal_ql)) + &
                                 kappa(3, 1)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, source)
            end if

            bcs_ht = localProps%auxiliar(1)                     ! downward flux at domain top

            ! solve radiative transfer equation along z
            call IR_RTE1_Global(localProps, nx*ny, nz, fdmi, source, b, flux_down, flux_up)
            ! call IR_RTE1_Local(localProps, nx*ny, nz, fdmi, p_source, p_b, p_flux_down, tmp2, p_flux_up)
            ! call IR_RTE1_Incremental(localProps, nx*ny, nz, fdmi, p_source, p_b, p_flux_down, p_flux_up)

            ! -----------------------------------------------------------------------
        case (TYPE_IR_BAND)
            ! t_ht(1:nx*ny) = s(isize_field - nx*ny + 1:isize_field, kmax, inb_scal_T)           ! save T at the top boundary

            p_flux_down = 0.0_wp                                        ! initialize for accumulation of conttributions from each band
            if (present(flux)) p_flux_up = 0.0_wp
            p_source = 0.0_wp
            do iband = 1, nbands
                b = sigma*s(:, inb_scal_T)**4.0_wp* &                   ! emission function, Stefan-Boltzmann law
                    (beta(1, iband) + s(:, inb_scal_T)*(beta(2, iband) + s(:, inb_scal_T)*beta(3, iband)))

                p_source(1:nx*ny*nz) = kappa(1, iband)*s(:, inb_scal_ql) + &               ! absorption coefficient
                                       kappa(2, iband)*(s(:, 2) - s(:, inb_scal_ql)) + &
                                       kappa(3, iband) ! calculate absorption coefficient into tmp_rad3
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, p_source) ! multiply by density
                end if

                bcs_ht(1:nx*ny) = localProps%auxiliar(iband)

                ! solve radiative transfer equation along z
                call IR_RTE1_Global(localProps, nx*ny, nz, fdmi, p_source, b, p_flux_down, p_flux_up)
                ! call IR_RTE1_Local(localProps, nx*ny, nz, fdmi, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), tmp2, p_b)
                ! call IR_RTE1_Incremental(localProps, nx*ny, nz, fdmi, tmp_rad(:, 3), tmp_rad(:, 5), tmp_rad(:, 4), p_b)
                if (present(flux)) then
                    flux_down = flux_down + p_flux_down
                    flux_up = flux_up + p_flux_up
                end if

                source = source + p_source

            end do

        end select

        ! if (present(flux)) then
        !     b = flux_up
        !     flux = flux_up - flux_down      ! net flux upwards
        ! end if

        nullify (p_tau)

    end subroutine Radiation_Infrared_Z

    !########################################################################
    ! Solve radiative transfer equation along 1 direction
    !########################################################################
    subroutine IR_RTE1_OnlyLiquid(localProps, nlines, ny, fdmi, a_source, flux_down, flux_up)
        type(radterm_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)      ! input as bulk absorption coefficent, output as source
        real(wp), intent(out), optional :: flux_down(nlines, ny), flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j

        ! #######################################################################
        ! calculate f_j = exp(-tau(z, zmax)/\mu)
        p_tau(:, ny) = 0.0_wp                                   ! boundary condition
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, a_source, p_tau, wrk2d)         ! recall this gives the negative of the integral
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), p_tau, BCS_MAX)
        ! call Int_Simpson_Biased_f(a_source, z%nodes(:), p_tau, BCS_MAX)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

        ! Calculate heating rate
        if (abs(localProps%auxiliar(2)) > 0.0_wp) then
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*(p_tau(:, j)*bcs_ht(1:nlines) &                     ! downward flux
                                                 + p_tau(:, 1)/p_tau(:, j)*bcs_hb(1:nlines))        ! upward flux
            end do
        else
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*p_tau(:, j)*bcs_ht(1:nlines)                        ! only downward
            end do
        end if

        ! Calculate flux, if necessary
        if (present(flux_up)) then
            do j = ny, 1, -1
                flux_down(:, j) = bcs_ht(1:nlines)*p_tau(:, j)                  ! downward flux
                flux_up(:, j) = bcs_hb(1:nlines)*p_tau(:, 1)/p_tau(:, j)        ! upward flux
            end do

        end if

        return
    end subroutine IR_RTE1_OnlyLiquid

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_Incremental(localProps, nlines, ny, fdmi, a_source, b, flux_down, flux_up)
        type(radterm_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)         ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)                ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: flux_down(nlines, ny)             ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(out), optional :: flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy
        real(wp), pointer :: p_wrk2d_1(:) => null()
        real(wp), pointer :: p_wrk2d_2(:) => null()

        ! #######################################################################
        p_wrk2d_1(1:nlines) => wrk2d(1:nlines, 1)
        p_wrk2d_2(1:nlines) => wrk2d(1:nlines, 2)

        ! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! transmission function I_{j-1,j}  = exp(-tau(z_{j-1}, z_j)/\mu)
        p_tau(:, 1) = 0.0_wp                                    ! boundary condition
        ! call FDM_Int1_Solve(nxz, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, p_tau, wrk2d)
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), p_tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), p_tau, BCS_MIN)
        do j = ny, 2, -1
            p_tau(:, j) = exp(p_tau(:, j - 1) - p_tau(:, j))
        end do
        p_tau(:, 1) = 1.0_wp        ! this value is not used

        ! ###################################################################
        ! downward flux; positive going down
        j = ny
        flux_down(1:nlines, j) = bcs_ht(1:nlines)

        do j = ny - 1, 1, -1
            ! Integral contribution from emission function using a trapezoidal rule
            p_wrk2d_2 = b(:, j + 1)
            p_wrk2d_1 = b(:, j)/p_tau(:, j + 1)
            p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(z%nodes(j + 1) - z%nodes(j))

            flux_down(:, j) = p_tau(:, j + 1)*(flux_down(:, j + 1) + p_wrk2d_1)
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and source
        j = 1
        a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)

        if (present(flux_up)) then                                          ! Save fluxes
            flux_up(:, j) = bcs_hb(1:nlines)                                ! upward flux

            do j = 2, ny
                ! Integral contribution from emission function using a trapezoidal rule
                p_wrk2d_1 = b(:, j - 1)
                p_wrk2d_2 = b(:, j)/p_tau(:, j)
                p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(z%nodes(j) - z%nodes(j - 1))

                bcs_hb = p_tau(:, j)*(bcs_hb(1:nlines) + p_wrk2d_1)
                a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)

                flux_up(:, j) = bcs_hb(1:nlines)
            end do

        else                                                                ! Do not save fluxes
            do j = 2, ny
                ! Integral contribution from emission function using a trapezoidal rule
                p_wrk2d_1 = b(:, j - 1)
                p_wrk2d_2 = b(:, j)/p_tau(:, j)
                p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(z%nodes(j) - z%nodes(j - 1))

                bcs_hb(1:nlines) = p_tau(:, j)*(bcs_hb(1:nlines) + p_wrk2d_1(1:nlines))
                a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)
            end do

        end if

        ! ###################################################################
        nullify (p_wrk2d_1, p_wrk2d_2)

        return
    end subroutine IR_RTE1_Incremental

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_Local(localProps, nlines, ny, fdmi, a_source, b, flux_down, tmp2, flux_up)
        type(radterm_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)         ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)                ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: flux_down(nlines, ny)             ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(inout) :: tmp2(nlines, ny)
        real(wp), intent(out), optional :: flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j, k
        real(wp) dummy
        real(wp), pointer :: p_flux_1(:) => null()
        real(wp), pointer :: p_flux_2(:) => null()

        ! #######################################################################
        p_flux_1(1:nlines) => wrk2d(1:nlines, 1)
        p_flux_2(1:nlines) => wrk2d(1:nlines, 2)

        ! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! transmission function I_{j-1,j} = exp(-tau(z_{j-1}, z_j)/\mu)
        p_tau(:, 1) = 0.0_wp                                    ! boundary condition
        ! call FDM_Int1_Solve(nxz, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, p_tau, wrk2d)
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), p_tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), p_tau, BCS_MIN)
        do j = ny, 2, -1
            p_tau(:, j) = exp(p_tau(:, j - 1) - p_tau(:, j))
        end do
        p_tau(:, 1) = 1.0_wp        ! this value is not used

        ! ###################################################################
        ! downward flux; positive going down
        j = ny
        p_flux_1(1:nlines) = bcs_ht(1:nlines)
        flux_down(:, j) = p_flux_1

        do j = ny - 1, 1, -1
            p_flux_1 = p_flux_1*p_tau(:, j + 1)                     ! accumulate transmission

            p_flux_2 = 1.0_wp                                       ! calculate emission; p_flux_2 is aux array in the next loop
            tmp2(:, j) = b(:, j)
            do k = j + 1, ny
                p_flux_2 = p_flux_2*p_tau(:, k)
                tmp2(:, k) = b(:, k)*p_flux_2
            end do
            call Int_Simpson_v(tmp2(:, j:ny), z%nodes(j:ny), flux_down(:, j))
            ! call Int_Trapezoidal_v(tmp2(:, j:ny), g%nodes(j:ny), flux_down(:, j))
            flux_down(:, j) = flux_down(:, j) + p_flux_1
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and net terms
        j = 1
        p_flux_1(1:nlines) = bcs_hb(1:nlines)
        a_source(:, j) = a_source(:, j)*(p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

        if (present(flux_up)) then                                  ! I need additional space to store fluxes
            flux_up(:, j) = p_flux_1                                ! upward flux

            do j = 2, ny
                p_flux_1 = p_flux_1*p_tau(:, j)

                p_flux_2 = 1.0_wp       ! used as auxiliary array in the next loop
                tmp2(:, j) = b(:, j)
                do k = j - 1, 1, -1
                    p_flux_2 = p_flux_2*p_tau(:, k + 1)
                    tmp2(:, k) = b(:, k)*p_flux_2
                end do
                call Int_Simpson_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                ! call Int_Trapezoidal_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                a_source(:, j) = a_source(:, j)*(p_flux_2 + p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

                flux_up(:, j) = p_flux_2 + p_flux_1                 ! upward flux
            end do

        else                ! we only need the source term
            do j = 2, ny
                p_flux_1 = p_flux_1*p_tau(:, j)                     ! accumulate transmission

                p_flux_2 = 1.0_wp                                   ! calculate emission; p_flux_2 is aux array in the next loop
                tmp2(:, j) = b(:, j)
                do k = j - 1, 1, -1
                    p_flux_2 = p_flux_2*p_tau(:, k + 1)
                    tmp2(:, k) = b(:, k)*p_flux_2
                end do
                call Int_Simpson_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                ! call Int_Trapezoidal_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                a_source(:, j) = a_source(:, j)*(p_flux_2 + p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

            end do

        end if

        ! #######################################################################
        nullify (p_flux_1, p_flux_2)

        return
    end subroutine IR_RTE1_Local

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_Global(localProps, nlines, ny, fdmi, a_source, b, flux_down, flux_up)
        type(radterm_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)     ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)            ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: flux_down(nlines, ny)    ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(inout) :: flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy

        ! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! ###################################################################
        ! downward flux; positive going down

        ! transmission function I_j = exp(-tau(z_j, zmax)/\mu)
        p_tau(:, ny) = 0.0_wp                                   ! boundary condition
        ! call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MIN)%rhs, a_source, p_tau, wrk2d)         ! recall this gives the negative of the integral
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), p_tau, BCS_MAX)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), p_tau, BCS_MAX)
        do j = ny, 1, -1
            p_tau(:, j) = exp(-p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

        flux_down = b/p_tau
        ! call Int_Trapezoidal_Increments_InPlace(flux_down, z%nodes(:), BCS_MAX)                   ! Calculate I_j = int_{x_{j}}^{x_{j+1}}
        call Int_Simpson_Biased_Increments_InPlace(flux_down, z%nodes(:), wrk2d(:, 1), BCS_MAX)   ! Calculate I_j = int_{x_{j}}^{x_{j+1}}
        j = ny
        wrk2d(1:nlines, 1) = 0.0_wp                                         ! accumulate emission; using wrk2d as aux array
        flux_down(:, j) = p_tau(:, j)*(bcs_ht(1:nlines) + wrk2d(1:nlines, 1))
        do j = ny - 1, 1, -1
            wrk2d(1:nlines, 1) = wrk2d(1:nlines, 1) + flux_down(:, j)            ! accumulate emission
            flux_down(:, j) = p_tau(:, j)*(bcs_ht(1:nlines) + wrk2d(1:nlines, 1))
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and net terms

        ! transmission function I_j = exp(-tau(zmin, z)/\mu)
        p_tau(:, 1) = 0.0_wp                                                ! boundary condition
        ! call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, p_tau, wrk2d)
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), p_tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), p_tau, BCS_MIN)
        do j = 1, ny
            p_tau(:, j) = exp(-p_tau(:, j))
        end do

        flux_up = b/p_tau
        ! call Int_Trapezoidal_Increments_InPlace(flux_up, z%nodes(:), BCS_MIN)                     ! Calculate I_j = int_{x_{j-1}}^{x_{j}}
        call Int_Simpson_Biased_Increments_InPlace(flux_up, z%nodes(:), wrk2d(:, 1), BCS_MIN)     ! Calculate I_j = int_{x_{j-1}}^{x_{j}}
        j = 1
        wrk2d(1:nlines, 1) = 0.0_wp                                         ! accumulate emission; using wrk2d as aux array
        flux_up(:, j) = p_tau(:, j)*(bcs_hb(1:nlines) + wrk2d(1:nlines, 1))
        !
        a_source(:, j) = a_source(:, j)*(flux_up(:, j) + flux_down(:, j)) - 2.0_wp*b(:, j)
        do j = 2, ny
            wrk2d(1:nlines, 1) = wrk2d(1:nlines, 1) + flux_up(:, j)         ! accumulate emission
            flux_up(:, j) = p_tau(:, j)*(bcs_hb(1:nlines) + wrk2d(1:nlines, 1))
            !
            a_source(:, j) = a_source(:, j)*(flux_up(:, j) + flux_down(:, j)) - 2.0_wp*b(:, j)
        end do

        return
    end subroutine IR_RTE1_Global

end module
