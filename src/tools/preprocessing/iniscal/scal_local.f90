#include "tlab_error.h"

module SCAL_LOCAL
    use TLab_Constants, only: wp, wi, pi_wp, big_wp, MAX_VARS
    use TLab_Constants, only: wfile, efile, lfile, tag_scal
    use TLab_Memory, only: imax, jmax, kmax, inb_scal, inb_txc
    use TLab_Time, only: itime
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: p_wrk2d, p_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_j
#endif
    use IO_Fields
    use TLab_Grid, only: x, y, z
    use Tlab_Background, only: sbg
    use Discrete, only: discrete_dt, Discrete_ReadBlock
    use Averages, only: AVG1V2D
    use Profiles, only: profiles_dt, Profiles_ReadBlock, Profiles_Calculate
    use Profiles, only: PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE, PROFILE_TANH_COS
    implicit none
    private

    public :: Iniscal_Initialize_Parameters
    public :: SCAL_FLUCTUATION_PLANE, SCAL_FLUCTUATION_VOLUME

    integer(wi), public :: flag_s                               ! Type of perturbation
    integer, parameter, public :: PERT_NONE = 0
    integer, parameter, public :: PERT_LAYER_BROADBAND = 1
    integer, parameter, public :: PERT_LAYER_DISCRETE = 2
    integer, parameter, public :: PERT_PLANE_BROADBAND = 4
    integer, parameter, public :: PERT_PLANE_DISCRETE = 5
    integer, parameter, public :: PERT_DELTA_BROADBAND = 6
    integer, parameter, public :: PERT_DELTA_DISCRETE = 7
    integer, parameter, public :: PERT_FLUX_BROADBAND = 8
    integer, parameter, public :: PERT_FLUX_DISCRETE = 9

    type(profiles_dt), public :: IniS(MAX_VARS)                 ! Geometry of perturbation of initial boundary condition

    real(wp), public :: norm_ini_radiation                      ! Scaling of perturbation
    integer(wi), public :: flag_mixture

    ! -------------------------------------------------------------------
    real(wp) :: norm_ini_s(MAX_VARS)                            ! Scaling of perturbation
    type(discrete_dt) :: fp                                     ! Discrete perturbation
    type(profiles_dt) :: prof_loc
    integer(wi) i, j, k
    integer(wi) im, idsp, jdsp
    real(wp) wx, wy, wx_1, wy_1

contains

    ! ###################################################################
    subroutine Iniscal_Initialize_Parameters(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        character(len=64) lstr
        integer(wi) idummy, is

        integer :: IniSvalid(5) = [PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE]

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'IniFields'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Scalar=<option>')
        call TLab_Write_ASCII(bakfile, '#NormalizeS=<values>')
        call TLab_Write_ASCII(bakfile, '#Mixture=<string>')

        call ScanFile_Char(bakfile, inifile, block, 'Scalar', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_s = 0
        else if (trim(adjustl(sRes)) == 'layerbroadband') then; flag_s = PERT_LAYER_BROADBAND
        else if (trim(adjustl(sRes)) == 'layerdiscrete') then; flag_s = PERT_LAYER_DISCRETE
        else if (trim(adjustl(sRes)) == 'planebroadband') then; flag_s = PERT_PLANE_BROADBAND
        else if (trim(adjustl(sRes)) == 'planediscrete') then; flag_s = PERT_PLANE_DISCRETE
        else if (trim(adjustl(sRes)) == 'deltabroadband') then; flag_s = PERT_DELTA_BROADBAND
        else if (trim(adjustl(sRes)) == 'deltadiscrete') then; flag_s = PERT_DELTA_DISCRETE
        else if (trim(adjustl(sRes)) == 'fluxbroadband') then; flag_s = PERT_FLUX_BROADBAND
        else if (trim(adjustl(sRes)) == 'fluxdiscrete') then; flag_s = PERT_FLUX_DISCRETE
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Scalar forcing type unknown')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call Profiles_ReadBlock(bakfile, inifile, block, 'IniS', IniS(1), 'gaussiansurface')
        if (trim(adjustl(sRes)) /= 'none') then
            IniS(2:) = IniS(1)
        else
            do is = 1, inb_scal
                write (lstr, *) is
                call Profiles_ReadBlock(bakfile, inifile, block, 'IniS'//trim(adjustl(lstr)), IniS(is), 'gaussiansurface')
            end do
        end if
        do is = 1, inb_scal
            if (.not. any(IniSvalid == IniS(is)%type)) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Undeveloped IniS type.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end do
        IniS(:)%delta = 1.0_wp
        IniS(:)%mean = 0.0_wp

        call ScanFile_Char(bakfile, inifile, block, 'NormalizeS', '-1.0', sRes)
        norm_ini_s(:) = 0.0_wp; idummy = inb_scal
        call LIST_REAL(sRes, idummy, norm_ini_s)
        if (idummy /= inb_scal) then            ! Consistency check
            if (idummy == 1) then
                norm_ini_s(2:) = norm_ini_s(1)
                call TLab_Write_ASCII(wfile, trim(adjustl(eStr))//'Using NormalizeS(1) for all scalars.')
            else
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'NormalizeS size does not match number of scalars.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        call ScanFile_Real(bakfile, inifile, block, 'NormalizeR', '0.0', norm_ini_radiation) ! Radiation field

        ! Additional parameters
        call ScanFile_Char(bakfile, inifile, block, 'Mixture', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_mixture = 0
        else if (trim(adjustl(sRes)) == 'equilibrium') then; flag_mixture = 1
        else if (trim(adjustl(sRes)) == 'loadfields') then; flag_mixture = 2
        end if

        ! ###################################################################
        block = 'Discrete'

        call Discrete_ReadBlock(bakfile, inifile, block, fp) ! Modulation type in fp%type
!   specific for this tool
        call ScanFile_Real(bakfile, inifile, block, 'Broadening', '-1.0', fp%parameters(1))
        call ScanFile_Real(bakfile, inifile, block, 'ThickStep', '-1.0', fp%parameters(2))

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        inb_txc = 0
        if (flag_s == PERT_LAYER_BROADBAND) inb_txc = max(inb_txc, 1)
        if (norm_ini_radiation /= 0.0_wp) inb_txc = max(inb_txc, 4)

        return
    end subroutine Iniscal_Initialize_Parameters

! ###################################################################
    subroutine SCAL_SHAPE(is, prof)
        integer(wi) is
        real(wp), intent(out) :: prof(:)

        ! -------------------------------------------------------------------
        real(wp) zr
        integer(wi) nz

#define zn(k) z%nodes(k)

        ! ###################################################################
        nz = size(prof)

        prof_loc = IniS(is)
        prof_loc%delta = 1.0_wp
        prof_loc%mean = 0.0_wp
        do k = 1, nz
            prof(k) = Profiles_Calculate(prof_loc, zn(k))
        end do

        select case (IniS(is)%type)
        case (PROFILE_GAUSSIAN_SURFACE) ! set perturbation and its normal derivative to zero at the boundaries
            do k = 1, nz
                zr = 0.5_wp*(zn(k) - zn(1))/IniS(is)%thick
                prof(k) = prof(k)*tanh(zr)**2
                zr = -0.5_wp*(zn(k) - zn(nz))/IniS(is)%thick
                prof(k) = prof(k)*tanh(zr)**2
            end do

        end select

#undef zn

        return
    end subroutine SCAL_SHAPE

! ###################################################################
    subroutine SCAL_FLUCTUATION_VOLUME(is, s, tmp)
        integer(wi) is
        real(wp), intent(out) :: s(imax, jmax, kmax)
        real(wp), intent(inout) :: tmp(imax, jmax, kmax)

        ! -------------------------------------------------------------------
        real(wp) dummy, amplify, params(0)

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; jdsp = ims_offset_j
#else
        idsp = 0; jdsp = 0
#endif

#define xn(i) x%nodes(i)
#define yn(j) y%nodes(j)

        call SCAL_SHAPE(is, wrk1d(1:kmax, 1))

        select case (flag_s)
        case (PERT_LAYER_BROADBAND)
            call IO_Read_Fields(trim(adjustl(tag_scal))//'rand', imax, jmax, kmax, itime, inb_scal, is, tmp, params)

            amplify = 0.0_wp
            do k = 1, kmax
                dummy = AVG1V2D(imax, jmax, kmax, k, 1, tmp)            ! Calculate mean
                p_wrk3d(:, :, k) = (tmp(:, :, k) - dummy)*wrk1d(k, 1)   ! Remove mean and apply shape function
            end do

        case (PERT_LAYER_DISCRETE)
            wx_1 = 2.0_wp*pi_wp/x%scale ! Fundamental wavelengths
            wy_1 = 2.0_wp*pi_wp/y%scale

            p_wrk2d = 0.0_wp
            do im = 1, fp%size
                wx = real(fp%modex(im), wp)*wx_1
                wy = real(fp%modey(im), wp)*wy_1
                do j = 1, jmax
                    p_wrk2d(:, j, 1) = p_wrk2d(:, j, 1) &
                                       + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*cos(wy*yn(jdsp + j) + fp%phasey(im))
                end do
            end do

            do k = 1, kmax
                p_wrk3d(:, :, k) = p_wrk2d(:, :, 1)*wrk1d(k, 1)
            end do

        end select

        if (norm_ini_s(is) > 0.0_wp) call SCAL_NORMALIZE(is, p_wrk3d)

        s = s + p_wrk3d

#undef xn
#undef yn

        return
    end subroutine SCAL_FLUCTUATION_VOLUME

    ! ###################################################################
    subroutine SCAL_FLUCTUATION_PLANE(is, s)
        integer(wi) is
        real(wp), intent(out) :: s(imax, jmax, kmax)

        ! -------------------------------------------------------------------
        real(wp) dummy
        real(wp) xcenter, ycenter, rcenter, amplify, params(0)

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; jdsp = ims_offset_j
#else
        idsp = 0; jdsp = 0
#endif

        ! ###################################################################
#define xn(i) x%nodes(i)
#define yn(j) z%nodes(k)
#define disp(i,k) p_wrk2d(i,k,1)

        disp(:, :) = 0.0_wp
        select case (flag_s)
        case (PERT_PLANE_BROADBAND, PERT_DELTA_BROADBAND, PERT_FLUX_BROADBAND)
            call IO_Read_Fields(trim(adjustl(tag_scal))//'rand', imax, jmax, 1, itime, 1, 1, disp(:, :), params)
            dummy = AVG1V2D(imax, jmax, 1, 1, 1, disp(:, :))      ! remove mean
            disp(:, :) = disp(:, :) - dummy

        case (PERT_PLANE_DISCRETE, PERT_DELTA_DISCRETE, PERT_FLUX_DISCRETE)
            wx_1 = 2.0_wp*pi_wp/x%scale ! Fundamental wavelengths
            wy_1 = 2.0_wp*pi_wp/y%scale

            do im = 1, fp%size
                wx = real(fp%modex(im), wp)*wx_1
                wy = real(fp%modey(im), wp)*wy_1

                if (fp%type == PROFILE_TANH_COS) then ! Smoothed step funtion Tanh(a*Cos(\xi/b))
                    if (fp%parameters(2) <= 0.0_wp) then
                        dummy = big_wp
                    else
                        dummy = 0.5_wp/(wx*fp%parameters(2))
                    end if
                    do j = 1, jmax
                        disp(:, j) = disp(:, j) &
                                     + fp%amplitude(im)*tanh(dummy*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*cos(wy*yn(jdsp + j) + fp%phasey(im)))
                    end do

                else
                    do j = 1, jmax
                        disp(:, j) = disp(:, j) &
                                     + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*cos(wy*yn(jdsp + j) + fp%phasey(im))
                    end do

                end if

            end do

        end select

        ! Modulation
        if (fp%type == PROFILE_GAUSSIAN .and. fp%parameters(1) > 0.0_wp) then
            do j = 1, jmax
                do i = 1, imax
                    xcenter = x%nodes(i + idsp) - x%scale*fp%phasex(1) - x%nodes(1)
                    if (y%size > 1) then
                        ycenter = y%nodes(j + jdsp) - y%scale*fp%phasey(1) - y%nodes(1)
                    else
                        ycenter = 0.0_wp
                    end if
                    rcenter = sqrt(xcenter**2 + ycenter**2)
                    amplify = exp(-0.5_wp*(rcenter/fp%parameters(1))**2)
                    disp(i, j) = disp(i, j)*amplify
                end do
            end do
        end if

        ! ###################################################################
        select case (flag_s)
        case (PERT_PLANE_BROADBAND, PERT_PLANE_DISCRETE)    ! Perturbation in the centerplane
            do j = 1, jmax
                do i = 1, imax

                    do k = 1, kmax
                        s(i, j, k) = Profiles_Calculate(sbg(is), z%nodes(k) - disp(i, j))
                    end do

                end do
            end do

        case (PERT_DELTA_BROADBAND, PERT_DELTA_DISCRETE)    ! Perturbation in the thickness
            prof_loc = sbg(is)

            do j = 1, jmax
                do i = 1, imax
                    prof_loc%thick = sbg(is)%thick + disp(i, j)

                    do k = 1, kmax
                        s(i, j, k) = Profiles_Calculate(prof_loc, z%nodes(k))
                    end do

                end do
            end do

        case (PERT_FLUX_BROADBAND, PERT_FLUX_DISCRETE)      ! Perturbation in the magnitude (constant derivative)
            prof_loc = sbg(is)

            do j = 1, jmax
                do i = 1, imax
                    prof_loc%delta = sbg(is)%delta + disp(i, j)
                    prof_loc%mean = (prof_loc%delta)*0.5_wp
                    if (sbg(is)%delta > 0) prof_loc%thick = prof_loc%delta/sbg(is)%delta*sbg(is)%thick

                    do k = 1, kmax
                        s(i, j, k) = Profiles_Calculate(prof_loc, z%nodes(k))
                    end do

                end do
            end do

        end select

#undef nx
#undef zn
#undef disp

        return
    end subroutine SCAL_FLUCTUATION_PLANE

    ! ###################################################################
    subroutine SCAL_NORMALIZE(is, s)
        integer(wi) is
        real(wp), intent(inout) :: s(imax, jmax, kmax)

        ! -------------------------------------------------------------------
        real(wp) dummy, amplify

        ! ###################################################################
        amplify = 0.0_wp                                ! Maximum across the layer
        do k = 1, kmax
            dummy = AVG1V2D(imax, jmax, kmax, k, 2, s)
            amplify = max(dummy, amplify)
        end do

        amplify = norm_ini_s(is)/sqrt(amplify)          ! Scaling factor to normalize to maximum rms

        s = s*amplify

        return
    end subroutine SCAL_NORMALIZE

end module SCAL_LOCAL
