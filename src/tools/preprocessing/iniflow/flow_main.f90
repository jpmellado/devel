program IniFlow
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, tag_flow
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Memory, only: inb_flow
    use TLab_Time, only: itime, rtime
    use TLab_Arrays
    use TLab_Pointers_3D, only: p_q
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use Thermodynamics, only: Thermo_Initialize
    ! use NavierStokes, only: nse_eqns, DNS_EQNS_COMPRESSIBLE, DNS_EQNS_TOTAL
    ! use Rotation, only: Rotation_Initialize
    use OPR_Partial, only: OPR_Partial_Initialize
    use TLab_Background, only: TLab_Initialize_Background, qbg
    use Profiles, only: Profiles_Calculate
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use OPR_Elliptic, only: OPR_Elliptic_Initialize
    use FLOW_LOCAL

    implicit none

    integer(wi) iq, k

    !########################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermo_Initialize(ifile)
    ! call Rotation_Initialize(ifile)

    call TLab_Consistency_Check()

    call Iniflow_Initialize_Parameters(ifile)

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call OPR_Partial_Initialize()
    if (flag_u /= 0) then
        call OPR_Fourier_Initialize()
        call OPR_Elliptic_Initialize(ifile)
        call OPR_Check()
    end if

    call TLab_Initialize_Background(ifile)
    if (IniK%relative) IniK%zmean = z%nodes(1) + z%scale*IniK%zmean_rel

    ! ###################################################################
    itime = 0; rtime = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing velocity fields.')

    ! Mean
    do iq = 1, 3
        do k = 1, kmax
            p_q(:, :, k, iq) = Profiles_Calculate(qbg(iq), z%nodes(k))
        end do
    end do
    ! if (coriolis%type == EQNS_COR_NORMALIZED) then
    !     ! call rotation()
    ! end if

    ! Fluctuation
    select case (flag_u)
    case (PERT_DISCRETE)
        call Iniflow_U_Discrete(txc(1, 1), txc(1, 2), txc(1, 3))
        q(1:isize_field, 1:3) = q(1:isize_field, 1:3) + txc(1:isize_field, 1:3)

    case (PERT_BROADBAND, PERT_BROADBAND_POTENTIAL, PERT_BROADBAND_VORTICITY)
        call Iniflow_U_Broadband(txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8))
        q(1:isize_field, 1:3) = q(1:isize_field, 1:3) + txc(1:isize_field, 1:3)

    end select

    ! ###################################################################
    ! Compressible formulation
    ! if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
    !     call TLab_Write_ASCII(lfile, 'Initializing pressure and density.')

    !     call PRESSURE_MEAN(p, T, s)
    !     call DENSITY_MEAN(rho, p, T, s, txc)

    !     if (flag_u /= 0) call PRESSURE_FLUCTUATION(q(1, 1), q(1, 2), q(1, 3), rho, p, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
    !     if (imixture /= 0) call IO_Read_Fields(trim(adjustl(tag_scal))//'ics', imax, jmax, kmax, itime, inb_scal, 0, s, params)
    !     if (flag_t /= 0) call DENSITY_FLUCTUATION(s, p, rho, txc(1, 1), txc(1, 2))

    !     ! Calculate specfic energy. Array s should contain the species fields at this point.
    !     call THERMO_THERMAL_TEMPERATURE(imax*jmax*kmax, s, p, rho, txc(:, 1))
    !     call THERMO_CALORIC_ENERGY(imax*jmax*kmax, s, txc(:, 1), e)

    ! end if

    ! ###################################################################
    io_header_q(1)%params(1) = rtime
    call IO_Write_Fields(trim(adjustl(tag_flow))//'ics', imax, jmax, kmax, itime, inb_flow, q, io_header_q)

    call TLab_Stop(0)

end program IniFlow
