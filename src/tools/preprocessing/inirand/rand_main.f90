program IniRand
    use TLab_Constants, only: wp, wi, tag_flow, tag_scal
    use TLab_Constants, only: ifile, gfile, lfile
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_Arrays, only: q, s, txc
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields, only: IO_Write_Fields, io_header_q, io_header_s
    use TLab_Grid
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use RAND_LOCAL

    implicit none

    ! -------------------------------------------------------------------
    integer(wi) iq, is

    ! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)

    call NavierStokes_Initialize_Parameters(ifile)

    call IniRand_Initialize_Parameters(ifile)

    ! ###################################################################
    call TLab_Initialize_Memory(__FILE__)

    call OPR_Fourier_Initialize()

    call OPR_Check()

    ! ###################################################################
    itime = 0; rtime = 0.0_wp

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Initializing random fields.')

    ! -------------------------------------------------------------------
    do iq = 1, inb_flow
        call RAND_FIELD(ucov(iq), q(:, iq), txc(:, 1), txc(:, 2), txc(:, 3))
    end do
    if (pdf%type == TYPE_DF_GAUSSIAN) then
        call RAND_COVARIANCE(ucov, q(:, 1), q(:, 2), q(:, 3))
    end if

    io_header_q(1)%params(1) = rtime
    call IO_Write_Fields(trim(adjustl(tag_flow))//'rand', imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))

    ! -------------------------------------------------------------------
    do is = 1, inb_scal
        call RAND_FIELD(ucov(is), s(:, is), txc(:, 1), txc(:, 2), txc(:, 3))
    end do

    io_header_s(:)%params(1) = rtime
    call IO_Write_Fields(trim(adjustl(tag_scal))//'rand', imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))

    call TLab_Stop(0)
end program IniRand
