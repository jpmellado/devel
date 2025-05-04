program VPARTIAL3D
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: gfile, ifile
    use TLab_Time, only: itime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: imax, jmax, kmax, inb_txc
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_err, ims_pro, ims_offset_i, ims_offset_j
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use FDM, only: g, fdm_dt, FDM_Initialize, FDM_CreatePlan
    use FDM_Derivative, only: FDM_COM4_DIRECT, FDM_COM6_JACOBIAN
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Grid
    use IO_Fields
    use OPR_Partial

    implicit none

    type(fdm_dt) :: g_loc
    real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht
    real(wp), dimension(:, :, :), pointer :: a, b, c, d, e, f
    real(wp), dimension(:, :, :), pointer :: u, du1_a, du2_a, du1_n, du2_n

    integer(wi) i, bcs(2, 2), idsp, jdsp
    integer(wi) type_of_problem
    real(wp) wk, x_0!, params(0) 

    ! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)

    inb_txc = 8

    call TLab_Initialize_Memory(__FILE__)

    allocate (bcs_ht(imax, kmax), bcs_hb(imax, kmax))

    u(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    du1_a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    du2_a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    du1_n(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    du2_n(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 7)

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    d(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    e(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    f(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 8)

    bcs = 0

    type_of_problem = 1     ! 1. order derivative
    ! type_of_problem = 2     ! 2. order derivative

    ! ###################################################################
    ! Define the function and analytic derivatives
    x_0 = 0.75_wp
    wk = 1.0_wp

#ifdef USE_MPI
    idsp = ims_offset_i; jdsp = ims_offset_j
#else
    idsp = 0; jdsp = 0
#endif

    do i = 1, imax
        ! single-mode
        u(i, :, :) = 1.0_wp + sin(2.0_wp*pi_wp/x%scale*wk*x%nodes(i + idsp)) ! + pi_wp/4.0_wp)
        du1_a(i, :, :) = (2.0_wp*pi_wp/x%scale*wk) &
                         *cos(2.0_wp*pi_wp/x%scale*wk*x%nodes(i + idsp))! + pi_wp/4.0_wp)
        du2_a(i, :, :) = -(2.0_wp*pi_wp/x%scale*wk)**2 &
                         *sin(2.0_wp*pi_wp/x%scale*wk*x%nodes(i + idsp))! + pi_wp/4.0_wp)
    end do

    ! do i = 1, jmax
    !     ! single-mode
    !     u(:, i, :) = 1.0_wp + sin(2.0_wp*pi_wp/y%scale*wk*y%nodes(i + jdsp)) ! + pi_wp/4.0_wp)
    !     du1_a(:, i, :) = (2.0_wp*pi_wp/y%scale*wk) &
    !                      *cos(2.0_wp*pi_wp/y%scale*wk*y%nodes(i + jdsp))! + pi_wp/4.0_wp)
    !     du2_a(:, i, :) = -(2.0_wp*pi_wp/y%scale*wk)**2 &
    !                      *sin(2.0_wp*pi_wp/y%scale*wk*y%nodes(i + jdsp))! + pi_wp/4.0_wp)
    ! end do

    ! do i = 1, kmax
    !     ! ! Gaussian
    !     u(:, :, i) = exp(-(z%nodes(i) - x_0*z%scale)**2/(2.0_wp*(z%scale/wk)**2))
    !     du1_a(:, :, i) = -(z%nodes(i) - x_0*z%scale)/(z%scale/wk)**2*u(:, :, i)
    !     du2_a(:, :, i) = -(z%nodes(i) - x_0*z%scale)/(z%scale/wk)**2*du1_a(:, :, i) &
    !                      - 1.0_wp/(z%scale/wk)**2*u(:, :, i)
    !     ! ! exponential
    !     ! u(:, i) = exp(-z%nodes(i)*wk)
    !     ! du1_a(:, i) = -wk*u(:, i)
    !     ! step
    !     ! u(:, i) = max(0.0_wp, (z%nodes(i) - z%nodes(kmax/2))*x_0)
    !     ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, z%nodes(i) - z%nodes(kmax/2)))*0.5_wp*x_0
    !     ! ! tanh
    !     ! u(:, i) = x_0*log(1.0_wp + exp((z%nodes(i) - z%nodes(kmax/2))/x_0))
    !     ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(z%nodes(i) - z%nodes(kmax/2))/x_0))
    !     ! ! Polynomial
    !     ! dummy = 4.0_wp
    !     ! u(:, i) = ((z%scale - z%nodes(i))/wk)**dummy
    !     ! du1_a(:, i) = -dummy*((z%scale - z%nodes(i))/wk)**(dummy - 1.0_wp)
    !     ! ! zero
    !     ! u(:, i) = 0.0_wp
    !     ! du1_a(:, i) = 0.0_wp
    !     ! ! delta-function
    !     ! u(:, i) = max(0.0_wp, 2.0_wp - real(i, wp))
    !     ! du1_a(:, i) = 0.0_wp
    !     ! du2_a(:, i) = 0.0_wp
    ! end do

    ! call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, f, params)

    ! ###################################################################
    select case (type_of_problem)
    case (1)
        call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, g(1), u, du2_n, du1_n)
        ! call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, g(2), u, du2_n, du1_n)
        ! call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, g(3), u, du2_n, du1_n)

        call check(du1_a, du1_n, txc(:, 1))

        call check(du2_a, du2_n, txc(:, 1))

        ! ###################################################################
    case (2)
        g_loc%der1%mode_fdm = FDM_COM6_JACOBIAN
        call FDM_CreatePlan(z, g_loc)
        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), f, c)
        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), c, a)
        call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, g_loc, f, a, c)
        call IO_Write_Fields('field.out1', imax, jmax, kmax, itime, 1, a)

        g_loc%der1%mode_fdm = FDM_COM4_DIRECT
        call FDM_CreatePlan(z, g_loc)
        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), f, d)
        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), d, b)
        call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, g_loc, f, b, d)
        call IO_Write_Fields('field.out2', imax, jmax, kmax, itime, 1, b)

        ! -------------------------------------------------------------------
        call check(a, b, txc(:, 1), 'field.dif')

        call check(c, d, txc(:, 1))

    end select

    call TLab_Stop(0)

    ! ###################################################################
contains
    subroutine check(a1, a2, dif, name)
        real(wp), intent(in) :: a1(*), a2(*)
        real(wp), intent(inout) :: dif(*)
        character(len=*), optional :: name

        real(wp) dummy, error
#ifdef USE_MPI
        real(wp) sum_mpi
#endif
        error = 0.0_wp
        dummy = 0.0_wp
        do i = 1, imax*jmax*kmax
            dif(i) = a2(i) - a1(i)
            error = error + dif(i)*dif(i)
            dummy = dummy + a1(i)*a1(i)
        end do
#ifdef USE_MPI
        sum_mpi = error
        call MPI_ALLREDUCE(sum_mpi, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        sum_mpi = dummy
        call MPI_ALLREDUCE(sum_mpi, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)

        if (ims_pro == 0) then
#endif
            write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
        end if
#endif

        if (present(name)) then
            call IO_Write_Fields(name, imax, jmax, kmax, itime, 1, dif)
        end if

        return
    end subroutine check

end program VPARTIAL3D
