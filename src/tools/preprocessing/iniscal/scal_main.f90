program IniScal
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, tag_scal
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal
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
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Background, only: TLab_Initialize_Background, sbg
    use Profiles, only: Profiles_Calculate
    use SCAL_LOCAL

    implicit none

    integer(wi) is, k

    ! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)

    call NavierStokes_Initialize_Parameters(ifile)

    call TLab_Consistency_Check()

    call IniScal_Initialize_Parameters(ifile)

    ! ###################################################################
    call TLab_Initialize_Memory(__FILE__)

    call TLab_Initialize_Background(ifile)
    do is = 1, size(IniS)
        if (IniS(is)%relative) IniS(is)%zmean = z%nodes(1) + z%scale*IniS(is)%zmean_rel
    end do

    ! ###################################################################
    itime = 0; rtime = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing scalar fields.')

    do is = 1, inb_scal
        ! Mean
        do k = 1, kmax
            p_s(:, :, k, is) = Profiles_Calculate(sbg(is), z%nodes(k))
        end do

        ! Fluctuation
        select case (flag_s)
        case (PERT_LAYER_BROADBAND, PERT_LAYER_DISCRETE)
            call SCAL_FLUCTUATION_VOLUME(is, s(:, is), txc)

        case (PERT_PLANE_BROADBAND, PERT_PLANE_DISCRETE, &
              PERT_DELTA_BROADBAND, PERT_DELTA_DISCRETE, PERT_FLUX_BROADBAND, PERT_FLUX_DISCRETE)
            call SCAL_FLUCTUATION_PLANE(is, s(:, is))

        end select

    end do

    ! ###################################################################
    io_header_s(:)%params(1) = rtime
    call IO_Write_Fields(trim(adjustl(tag_scal))//'ics', imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))

    call TLab_Stop(0)

end program IniScal
