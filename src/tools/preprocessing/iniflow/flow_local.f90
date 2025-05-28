#include "tlab_error.h"

module FLOW_LOCAL
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: efile, lfile, wfile, tag_flow
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Memory, only: inb_wrk2d, inb_txc
    use TLab_Time, only: itime!, rtime
    use TLab_Pointers_3D, only: p_wrk1d, p_wrk2d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_j
#endif
    use IO_Fields
    use TLab_Grid, only: x, y, z
    use FDM, only: g
    use Averages, only: AVG1V2D
    use Profiles, only: profiles_dt, Profiles_ReadBlock, Profiles_Calculate
    use Profiles, only: PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE, PROFILE_PARABOLIC_SURFACE
    use OPR_Partial
    use OPR_Elliptic
    use Discrete, only: discrete_dt, Discrete_ReadBlock
    use FI_VECTORCALCULUS
    implicit none
    private

    public :: Iniflow_Initialize_Parameters
    public :: Iniflow_U_Broadband, Iniflow_U_Discrete
    ! public :: DENSITY_FLUCTUATION, PRESSURE_FLUCTUATION ! Only in compressible formulation

    integer(wi), public :: flag_u               ! Type of perturbation in velocity
    integer, parameter, public :: PERT_NONE = 0
    integer, parameter, public :: PERT_DISCRETE = 1
    integer, parameter, public :: PERT_BROADBAND = 2
    integer, parameter, public :: PERT_BROADBAND_VORTICITY = 3
    integer, parameter, public :: PERT_BROADBAND_POTENTIAL = 4

    type(profiles_dt), public :: IniK           ! Geometry of perturbation of initial boundary condition

    ! integer(wi), public :: flag_t               ! Type of perturbation in thermodynamic fields

    ! -------------------------------------------------------------------
    integer(wi) :: ibc_pert                     ! BCs at kmin/kmax: 0, No-Slip/No-Slip
    !                                                               1, Free-Slip/No-Slip
    !                                                               2, No-Slip/Free-Slip
    !                                                               3, Free-Slip/Free-Slip
    logical :: RemoveDilatation

    real(wp) :: norm_ini_u                      ! Scaling of perturbation
    ! real(wp) :: norm_ini_p
    type(discrete_dt) :: fp                     ! Discrete perturbation
    integer(wi) im, idsp, jdsp
    real(wp) wx, wy, wx_1, wy_1

contains

    ! ###################################################################
    subroutine Iniflow_Initialize_Parameters(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer :: IniKvalid(6) = [PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, &
                                   PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE, PROFILE_PARABOLIC_SURFACE]

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Inifields'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Velocity=<VelocityDiscrete/VelocityBroadband/PotentialBroadband/VorticityBroadband>')
        call TLab_Write_ASCII(bakfile, '#Temperature=<option>')
        call TLab_Write_ASCII(bakfile, '#ForceDilatation=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#NormalizeK=<value>')
        call TLab_Write_ASCII(bakfile, '#NormalizeP=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Velocity', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_u = PERT_NONE
        else if (trim(adjustl(sRes)) == 'velocitydiscrete') then; flag_u = PERT_DISCRETE
        else if (trim(adjustl(sRes)) == 'velocitybroadband') then; flag_u = PERT_BROADBAND
        else if (trim(adjustl(sRes)) == 'vorticitybroadband') then; flag_u = PERT_BROADBAND_VORTICITY
        else if (trim(adjustl(sRes)) == 'potentialbroadband') then; flag_u = PERT_BROADBAND_POTENTIAL
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Velocity forcing type unknown')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        RemoveDilatation = .true.
        call ScanFile_Char(bakfile, inifile, block, 'ForceDilatation', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'no') RemoveDilatation = .false.

        call Profiles_ReadBlock(bakfile, inifile, block, 'IniK', IniK)
        if (.not. any(IniKvalid == IniK%type)) then
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Undeveloped IniK type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if
        IniK%delta = 1.0_wp
        IniK%mean = 0.0_wp

        call ScanFile_Real(bakfile, inifile, block, 'NormalizeK', '-1.0', norm_ini_u)

        ! ! Compressible formulation
        ! call ScanFile_Char(bakfile, inifile, block, 'Temperature', 'None', sRes)
        ! if (trim(adjustl(sRes)) == 'none') then; flag_t = 0
        ! else if (trim(adjustl(sRes)) == 'planediscrete') then; flag_t = PERT_DISCRETE
        ! else if (trim(adjustl(sRes)) == 'planebroadband') then; flag_t = PERT_BROADBAND
        ! else
        !     call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Temperature forcing type unknown')
        !     call TLab_Stop(DNS_ERROR_OPTION)
        ! end if

        ! call ScanFile_Real(bakfile, inifile, block, 'NormalizeP', '-1.0', norm_ini_p)

        ! ###################################################################
        block = 'BoundaryConditions'      ! Boundary conditions for the perturbation

        ibc_pert = 0
        call ScanFile_Char(bakfile, inifile, block, 'VelocityKmin', 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then
        else if (trim(adjustl(sRes)) == 'noslip') then
        else if (trim(adjustl(sRes)) == 'freeslip') then; ibc_pert = ibc_pert + 1
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'VelocityKmin.')
            call TLab_Stop(DNS_ERROR_IBC)
        end if
        call ScanFile_Char(bakfile, inifile, block, 'VelocityKmax', 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then
        else if (trim(adjustl(sRes)) == 'noslip') then
        else if (trim(adjustl(sRes)) == 'freeslip') then; ibc_pert = ibc_pert + 2
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'VelocityKmax.')
            call TLab_Stop(DNS_ERROR_IBC)
        end if

        ! ###################################################################
        block = 'Discrete'

        call Discrete_ReadBlock(bakfile, inifile, block, fp) ! Modulation type in fp%type
        ! specific for this tool
        call ScanFile_Real(bakfile, inifile, block, 'Broadening', '-1.0', fp%parameters(1))
        call ScanFile_Real(bakfile, inifile, block, 'ThickStep', '-1.0', fp%parameters(2))

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        inb_wrk2d = max(inb_wrk2d, 3)

        inb_txc = 2
        if (flag_u /= PERT_NONE) inb_txc = max(inb_txc, 8)

        return
    end subroutine Iniflow_Initialize_Parameters

    ! ###################################################################
    ! ###################################################################
    subroutine Iniflow_U_Discrete(u, v, w)
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
        use Thermo_Anelastic, only: ribackground, Thermo_Anelastic_Weight_InPlace
        real(wp), dimension(imax, jmax, kmax), intent(out) :: u, v, w

        ! -------------------------------------------------------------------
        real(wp) factorx, factory
        integer(wi) j, k

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; jdsp = ims_offset_j
#else
        idsp = 0; jdsp = 0
#endif

#define xn(i) x%nodes(i)
#define yn(j) y%nodes(j)

        call FLOW_SHAPE(p_wrk1d)

        wx_1 = 2.0_wp*pi_wp/x%scale ! Fundamental wavelengths
        wy_1 = 2.0_wp*pi_wp/y%scale

        p_wrk2d = 0.0_wp
        do im = 1, fp%size
            wx = real(fp%modex(im), wp)*wx_1
            wy = real(fp%modey(im), wp)*wy_1

            ! Factor to impose solenoidal constraint
            if (fp%modex(im) == 0 .and. fp%modey(im) == 0) then
                exit
            elseif (fp%modey(im) == 0) then
                factorx = 1.0_wp/wx; factory = 0.0_wp
            elseif (fp%modex(im) == 0) then
                factorx = 0.0_wp; factory = 1.0_wp/wy
            else
                factorx = 0.5_wp/wx; factory = 0.5_wp/wy
            end if

            do j = 1, jmax
                p_wrk2d(:, j, 3) = p_wrk2d(:, j, 3) + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                   *cos(wy*yn(jdsp + j) + fp%phasey(im))
                p_wrk2d(:, j, 1) = p_wrk2d(:, j, 1) + fp%amplitude(im)*sin(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                   *cos(wy*yn(jdsp + j) + fp%phasey(im))*factorx
                p_wrk2d(:, j, 2) = p_wrk2d(:, j, 2) + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                   *sin(wy*yn(jdsp + j) + fp%phasey(im))*factory
            end do

        end do

        do k = 1, kmax
            u(:, :, k) = p_wrk2d(:, :, 1)*p_wrk1d(k, 2)
            v(:, :, k) = p_wrk2d(:, :, 2)*p_wrk1d(k, 2)
            w(:, :, k) = p_wrk2d(:, :, 3)*p_wrk1d(k, 1)
        end do

        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, u)
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, v)
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, w)
        end if

        if (norm_ini_u >= 0.0_wp) call FLOW_NORMALIZE(u, v, w)

#undef xn
#undef yn

        return
    end subroutine Iniflow_U_Discrete

    ! ###################################################################
    subroutine Iniflow_U_Broadband(u, v, w, ax, ay, az, tmp4, tmp5)
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
        use Thermo_Anelastic, only: rbackground, ribackground, Thermo_Anelastic_Weight_InPlace
        use FI_VECTORCALCULUS

        real(wp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w
        real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: ax, ay, az, tmp4, tmp5

        ! -------------------------------------------------------------------
        integer ibc_loc
        real(wp) dummy, params(0)
        integer(wi) k

        real(wp), allocatable :: bcs_hb(:), bcs_ht(:)

        ! ###################################################################
        call IO_Read_Fields(trim(adjustl(tag_flow))//'rand', imax, jmax, kmax, itime, 3, 1, u, params)
        call IO_Read_Fields(trim(adjustl(tag_flow))//'rand', imax, jmax, kmax, itime, 3, 2, v, params)
        call IO_Read_Fields(trim(adjustl(tag_flow))//'rand', imax, jmax, kmax, itime, 3, 3, w, params)

        do k = 1, kmax   ! Remove mean
            dummy = AVG1V2D(imax, jmax, kmax, k, 1, u)
            u(:, :, k) = u(:, :, k) - dummy
            dummy = AVG1V2D(imax, jmax, kmax, k, 1, v)
            v(:, :, k) = v(:, :, k) - dummy
            dummy = AVG1V2D(imax, jmax, kmax, k, 1, w)
            w(:, :, k) = w(:, :, k) - dummy
        end do

        call FLOW_SHAPE(p_wrk1d)

        ! ###################################################################
        select case (flag_u)
        case (PERT_BROADBAND)                           ! Velocity u given
            do k = 1, kmax
                u(:, :, k) = u(:, :, k)*p_wrk1d(k, 2)
                v(:, :, k) = v(:, :, k)*p_wrk1d(k, 2)
                w(:, :, k) = w(:, :, k)*p_wrk1d(k, 1)
            end do

        case (PERT_BROADBAND_POTENTIAL)                 ! Velocity potential u given, calculate u = rot(u)
            do k = 1, kmax
                ax(:, :, k) = u(:, :, k)*p_wrk1d(k, 1)  ! Horizontal components of vector potential give vertical velocity
                ay(:, :, k) = v(:, :, k)*p_wrk1d(k, 1)
                az(:, :, k) = w(:, :, k)*p_wrk1d(k, 2)
            end do

            ! Cannot use fi_curl. I need to impose BCs to zero to get zero velocity there
            ibc_loc = 0
            if (any([BCS_DD, BCS_DN] == ibc_pert)) ibc_loc = ibc_loc + 1 ! no-slip at ymin
            if (any([BCS_DD, BCS_ND] == ibc_pert)) ibc_loc = ibc_loc + 2 ! no-slip at ymax
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), az, u)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), ay, tmp4, ibc=ibc_loc)
            u = u - tmp4
            if (y%size > 1) then
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), ax, v, ibc=ibc_loc)
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), az, tmp4)
            end if
            v = v - tmp4
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), ay, w)
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), ax, tmp4)
            w = w - tmp4

        case (PERT_BROADBAND_VORTICITY)                 ! Vorticity given, solve lap(u) = - rot(vort), vort = rot(u)
            if (allocated(bcs_hb)) deallocate (bcs_hb)
            if (allocated(bcs_ht)) deallocate (bcs_ht)
            allocate (bcs_hb(imax*jmax), bcs_ht(imax*jmax))

            call FI_CURL(imax, jmax, kmax, u, v, w, ax, ay, az, tmp4)
            do k = 1, kmax
                ax(:, :, k) = -ax(:, :, k)*p_wrk1d(k, 2)
                ay(:, :, k) = -ay(:, :, k)*p_wrk1d(k, 2)
                az(:, :, k) = -az(:, :, k)*p_wrk1d(k, 1)
            end do
            call FI_CURL(imax, jmax, kmax, ax, ay, az, u, v, w, tmp4)

            ! Solve lap(u) = - (rot(vort))_x
            bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
            call OPR_Poisson(imax, jmax, kmax, ibc_pert, u, tmp4, tmp5, bcs_hb, bcs_ht)

            ! Solve lap(v) = - (rot(vort))_y with no penetration bcs
            if (y%size > 1) then
                bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
                call OPR_Poisson(imax, jmax, kmax, ibc_pert, v, tmp4, tmp5, bcs_hb, bcs_ht)
            end if

            ! Solve lap(w) = - (rot(vort))_z
            bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
            call OPR_Poisson(imax, jmax, kmax, BCS_DD, w, tmp4, tmp5, bcs_hb, bcs_ht)

        end select

        ! ###################################################################
        if (RemoveDilatation) then              ! Remove dilatation (vort might not be a vorticity field because it was not solenoidal)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, u)
                call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, v)
                call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, w)
            end if
            call FI_SOLENOIDAL(imax, jmax, kmax, u, v, w, ax, ay, az)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, u)
                call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, v)
                call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, w)
            end if
        end if

        if (y%size == 1) v = 0.0_wp          ! Impose zero spanwise velocity in 2D case

        if (norm_ini_u >= 0.0_wp) call FLOW_NORMALIZE(u, v, w)

        return
    end subroutine Iniflow_U_Broadband

    ! ###################################################################
    subroutine FLOW_SHAPE(profs)
        real(wp), intent(inout) :: profs(kmax, 2)

        ! -------------------------------------------------------------------
        real(wp) zr
        integer(wi) k

#define zn(k) z%nodes(k)

        ! ###################################################################
        do k = 1, kmax                                              ! Wall-normal velocity
            profs(k, 1) = Profiles_Calculate(IniK, zn(k))
        end do
        call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), profs(:, 1), profs(:, 2))
        profs(:, 2) = -profs(:, 2)                                  ! Negative of the derivative of f, wall-parallel velocity

        select case (IniK%type)
        case (PROFILE_PARABOLIC_SURFACE)
            ! Zero wall-parallel velocity for no-slip condition, multiply by parabolic again, f=f*f
            profs(:, 2) = 2.0_wp*profs(:, 2)*profs(:, 1)            ! Wall-parallel velocity
            profs(:, 1) = profs(:, 1)**2.0_wp                       ! Wall-normal velocity

        case (PROFILE_GAUSSIAN_SURFACE)
            ! Zero wall-normal derivative of wall-parallel velocity for free-slip and potentialvelocity mode, f=f*tanh
            if (any([BCS_DD, BCS_DN] == ibc_pert)) then            ! no-slip at kmin
                do k = 1, kmax
                    zr = 0.5_wp*(zn(k) - zn(1))/IniK%thick
                    profs(k, 2) = profs(k, 2)*tanh(zr)**2 - &       ! Wall-parallel velocity
                                  profs(k, 1)*tanh(zr)/cosh(zr)**2/IniK%thick
                    profs(k, 1) = profs(k, 1)*tanh(zr)**2           ! Wall-normal velocity
                end do
            end if

            if (any([BCS_DD, BCS_ND] == ibc_pert)) then            ! no-slip at kmax
                do k = 1, kmax
                    zr = 0.5_wp*(zn(kmax) - zn(k))/IniK%thick
                    profs(k, 2) = profs(k, 2)*tanh(zr)**2 + &       ! Wall-parallel velocity
                                  profs(k, 1)*tanh(zr)/cosh(zr)**2/IniK%thick
                    profs(k, 1) = profs(k, 1)*tanh(zr)**2           ! Wall-normal velocity
                end do
            end if

        end select

#undef zn

        return
    end subroutine FLOW_SHAPE

    ! ###################################################################
    subroutine FLOW_NORMALIZE(u, v, w)
        real(wp), dimension(imax, jmax, kmax) :: u, v, w

        ! -------------------------------------------------------------------
        real(wp) dummy, amplify
        integer(wi) k

        ! ###################################################################
        amplify = 0.0_wp                                ! Maximum across the layer
        do k = 1, kmax
            dummy = AVG1V2D(imax, jmax, kmax, k, 2, u) + &
                    AVG1V2D(imax, jmax, kmax, k, 2, v) + &
                    AVG1V2D(imax, jmax, kmax, k, 2, w)
            amplify = max(dummy, amplify)
        end do
        amplify = 0.5_wp*amplify

        amplify = sqrt(norm_ini_u/amplify)              ! Scaling factor to normalize to maximum TKE

        u = u*amplify
        v = v*amplify
        w = w*amplify

        return
    end subroutine FLOW_NORMALIZE

!     ! ###################################################################
!     !# Perturbation of the thermodynamic fields by a displacement of the reference center plane.
!     !# Array s enters with the scalar total field, including fluctuations.
!     !# Only used in compressible formulation
!     !# Together discrete and broadband in one procedure
!     !########################################################################
!     subroutine DENSITY_FLUCTUATION(s, p, rho, T, h)
!         use Thermodynamics, only: imixture, MIXT_TYPE_AIRWATER

!         real(wp), dimension(imax, jmax, kmax) :: T, h, rho, p
!         real(wp), dimension(imax, jmax, kmax, *) :: s

!         ! -------------------------------------------------------------------
!         real(wp) dummy
!         real(wp) xcenter, amplify, params(0)
!         type(profiles_dt) :: prof_loc

!         ! ###################################################################

! #ifdef USE_MPI
!         idsp = ims_offset_i; jdsp = ims_offset_j
! #else
!         idsp = 0; jdsp = 0
! #endif

!         ! ###################################################################
!         ! Center plane displacement
!         ! ###################################################################
! #define yn(j) y%nodes(j)
! #define xn(i) x%nodes(i)
! #define zn(k) z%nodes(k)
! #define disp(i,k) p_wrk2d(i,k,1)

!         disp(:, :) = 0.0_wp

!         select case (flag_t)
!         case (PERT_BROADBAND)
!             call IO_Read_Fields('scal.rand', imax, 1, kmax, itime, 1, 1, disp(:, :), params)
!             dummy = AVG1V2D(imax, 1, kmax, 1, 1, disp(:, :))     ! remove mean
!             disp(:, :) = disp(:, :) - dummy

!         case (PERT_DISCRETE)
!             wx_1 = 2.0_wp*pi_wp/x%scale ! Fundamental wavelengths
!             wy_1 = 2.0_wp*pi_wp/y%scale

!             do im = 1, fp%size
!                 wx = real(fp%modex(im), wp)*wx_1
!                 wy = real(fp%modey(im), wp)*wy_1

!                 do k = 1, kmax
!                     disp(:, k) = disp(:, k) + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
!                                  *cos(wy*yn(jdsp + k) + fp%phasey(im))
!                 end do

!             end do

!         end select

!         ! -------------------------------------------------------------------
!         ! Modulation
!         ! -------------------------------------------------------------------
!         if (fp%type == PROFILE_GAUSSIAN .and. fp%parameters(1) > 0.0_wp) then
!             do k = 1, kmax
!                 do i = 1, imax
!                     xcenter = xn(i) - x%scale*0.5_wp - xn(1)
!                     amplify = exp(-0.5_wp*(xcenter/fp%parameters(1))**2)
!                     disp(i, k) = disp(i, k)*amplify
!                 end do
!             end do
!         end if

!         ! ###################################################################
!         ! Perturbation in the thermodynamic fields
!         ! ###################################################################
!         if (tbg%type /= PROFILE_NONE) then
!             prof_loc = tbg
!             do k = 1, kmax
!                 do i = 1, imax
!                     prof_loc%ymean = tbg%ymean + disp(i, k)
!                     prof_loc%delta = tbg%delta + (tbg%uslope - tbg%lslope)*disp(i, k)*y%scale
!                     prof_loc%mean = tbg%mean + 0.5_wp*(tbg%uslope + tbg%lslope)*disp(i, k)*y%scale
!                     do j = 1, jmax
!                         T(i, j, k) = Profiles_Calculate(prof_loc, yn(j))
!                     end do
!                 end do
!             end do

!             if (imixture == MIXT_TYPE_AIRWATER) then
!                 call THERMO_AIRWATER_PT(imax*jmax*kmax, s, p, T)
!             end if

!             call THERMO_THERMAL_DENSITY(imax*jmax*kmax, s, p, T, rho)

!         end if

!         if (hbg%type /= PROFILE_NONE) then
!             prof_loc = hbg
!             do k = 1, kmax
!                 do i = 1, imax
!                     prof_loc%ymean = hbg%ymean + disp(i, k)
!                     prof_loc%delta = hbg%delta + (hbg%uslope - hbg%lslope)*disp(i, k)*y%scale
!                     prof_loc%mean = hbg%mean + 0.5_wp*(hbg%uslope + hbg%lslope)*disp(i, k)*y%scale
!                     do j = 1, jmax
!                         h(i, j, k) = Profiles_Calculate(prof_loc, yn(j))
!                     end do
!                 end do
!             end do

!             if (imixture == MIXT_TYPE_AIRWATER) then
!                 call THERMO_AIRWATER_PH_RE(imax*jmax*kmax, s, p, h, T)
!             end if

!             call THERMO_THERMAL_DENSITY(imax*jmax*kmax, s, p, T, rho)

!         end if

! #undef xn
! #undef yn
! #undef zn
! #undef disp

!         return
!     end subroutine DENSITY_FLUCTUATION

!     ! ###################################################################
!     !# solve Poisson equation for p', nabla^2 p' = d/dx_i d/dx_j (rho_0 u_i u_j),
!     !# assuming p/rho^\gamma0 constant(Homentropic conditions)
!     ! ###################################################################
!     subroutine PRESSURE_FLUCTUATION(u, v, w, rho, p, pprime, txc1, txc2, txc3, txc4)
!         use Thermodynamics, only: gamma0

!         real(wp), dimension(imax, jmax, kmax), intent(in) :: u, v, w
!         real(wp), dimension(imax, jmax, kmax), intent(inout) :: rho, p, pprime
!         real(wp), dimension(imax, jmax, kmax), intent(inout) :: txc1, txc2, txc3, txc4

!         ! -------------------------------------------------------------------
!         integer(wi) bcs(2, 2)

!         real(wp), allocatable :: bcs_hb(:), bcs_ht(:)

!         ! ###################################################################
!         ! Calculate RHS d/dx_i d/dx_j (u_i u_j), stored in txc4

!         ! terms with u
!         txc1 = rho*u*u; txc2 = rho*u*v; txc3 = rho*u*w
!         call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc3, txc4)
!         call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc2, txc3)
!         call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc1, txc2)
!         txc2 = 2.0_wp*(txc4 + txc3) + txc2

!         call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc2, txc4)

!         ! terms with v
!         txc1 = rho*v*v; txc2 = rho*v*w
!         call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc2, txc3)
!         call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc1, txc2)
!         txc2 = txc2 + 2.0_wp*txc3

!         call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc2, txc1)
!         txc4 = txc4 + txc1

!         ! terms with w
!         txc1 = rho*w*w
!         call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc1, txc2)

!         call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc2, txc1)
!         txc4 = txc4 + txc1

!         ! Solve Poisson equation; pprime contains fluctuating p' (BCs are equal to zero!)
!         if (g(1)%periodic .and. g(3)%periodic) then ! Doubly periodic in xOz
!             pprime = -txc4          ! change of forcing term sign

!             allocate (bcs_hb(imax*kmax), bcs_ht(imax*kmax))
!             bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
!             call OPR_Poisson(imax, jmax, kmax, 0, pprime, txc1, txc2, bcs_hb, bcs_ht)
!         else                                      ! General treatment
!             ! Undevelop
!         end if

!         ! An amplification factor norm_ini_p is allowed as in previous versions
!         rho = (norm_ini_p*pprime/p/gamma0 + 1.0_wp)*rho  ! isentropic relation p'/p = \gamma \rho'/rho
!         p = norm_ini_p*pprime + p

!         return
!     end subroutine PRESSURE_FLUCTUATION

end module FLOW_LOCAL
