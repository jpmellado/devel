#include "tlab_error.h"

#define C_FILE_LOC "INISCAL"

program INISCAL

    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, tag_scal
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal !inb_scal_array, inb_txc, isize_field, inb_wrk2d
    use TLab_Arrays
    use TLab_Pointers_3D, only: p_s
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    ! use FDM, only: fdm_Int0
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    ! use Thermodynamics
    use TLab_Background, only: TLab_Initialize_Background, sbg
    ! use Gravity, only: Gravity_Initialize
    ! use Rotation, only: Rotation_Initialize
    ! use Radiation, only: Radiation_Initialize, infraredProps, Radiation_Infrared_Y, radterm_dt, TYPE_RAD_NONE
    ! use LargeScaleForcing, only: LargeScaleForcing_Initialize
    ! use THERMO_AIRWATER
    ! use Thermo_Anelastic
    use Profiles, only: Profiles_Calculate
    use SCAL_LOCAL

    implicit none

    ! -------------------------------------------------------------------
    integer(wi) is, k
    ! type(radterm_dt) localInfraredProps

    ! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    ! call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    ! call Thermodynamics_Initialize_Parameters(ifile)
    ! call Radiation_Initialize(ifile)

    call TLab_Consistency_Check()

    call Iniscal_Initialize_Parameters(ifile)

    ! ###################################################################
    call TLab_Initialize_Memory(C_FILE_LOC)

    call TLab_Initialize_Background(ifile)
    do is = 1, size(IniS)
        if (IniS(is)%relative) IniS(is)%zmean = z%nodes(1) + z%scale*IniS(is)%zmean_rel
    end do

    ! ###################################################################
    itime = 0; rtime = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing scalars.')

    do is = 1, inb_scal
        do k = 1, kmax
            p_s(:, :, k, is) = Profiles_Calculate(sbg(is), z%nodes(k))
        end do

        select case (flag_s)
        case (PERT_LAYER_BROADBAND, PERT_LAYER_DISCRETE)
            call SCAL_FLUCTUATION_VOLUME(is, s(:, is), txc)

        case (PERT_PLANE_BROADBAND, PERT_PLANE_DISCRETE, &
              PERT_DELTA_BROADBAND, PERT_DELTA_DISCRETE, PERT_FLUX_BROADBAND, PERT_FLUX_DISCRETE)
            call SCAL_FLUCTUATION_PLANE(is, s(:, is))

        end select

    end do

    ! ! ###################################################################
    ! ! Initial liquid in equilibrium; overwrite previous values
    ! if (flag_mixture == 1 .or. inb_scal_array > inb_scal) then
    !     select case (imode_thermo)
    !     case (THERMO_TYPE_ANELASTIC)
    !         if (imixture == MIXT_TYPE_AIRWATER) then
    !             call Thermo_Anelastic_PH(imax, jmax, kmax, s(1, 2), s(1, 1))
    !         end if

    !     case (THERMO_TYPE_LINEAR)
    !         if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
    !             call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, inb_scal_array))
    !         end if

    !     end select

    ! end if

    ! ! ###################################################################
    ! ! Initial radiation effect as an accumulation during a certain interval of time
    ! if (infraredProps%type /= TYPE_RAD_NONE .and. norm_ini_radiation /= 0.0_wp) then
    !     norm_ini_radiation = norm_ini_radiation/infraredProps%auxiliar(1)
    !     localInfraredProps = infraredProps
    !     localInfraredProps%auxiliar(:) = localInfraredProps%auxiliar(:)*norm_ini_radiation

    !     do is = 1, inb_scal
    !         if (localInfraredProps%active(is)) then
    !             call Radiation_Infrared_Y(localInfraredProps, imax, jmax, kmax, fdm_Int0, s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
    !             s(1:isize_field, is) = s(1:isize_field, is) + txc(1:isize_field, 1)
    !         end if
    !     end do

    ! end if

    ! ###################################################################
    io_header_s(:)%params(1) = rtime
    call IO_Write_Fields(trim(adjustl(tag_scal))//'ics', imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))

    call TLab_Stop(0)
end program INISCAL
