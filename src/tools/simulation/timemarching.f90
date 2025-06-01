#include "tlab_error.h"

!########################################################################
!#
!# Runge-Kutta explicit 3th order from Williamson 1980
!# Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
!# Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
!#
!########################################################################
module TimeMarching

    use TLab_Constants, only: efile, wfile, wp, wi, big_wp
    use TLab_WorkFlow, only: flow_on, scal_on
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Memory, only: isize_wrk1d, isize_wrk2d, isize_wrk3d
    use TLab_Memory, only: isize_txc_field
    use TLab_Time, only: rtime
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    ! use PARTICLE_VARS
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS
#endif
    use TLab_Grid, only: x, y, z
    use NavierStokes, only: nse_eqns, DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC
    use NavierStokes, only: nse_advection, EQNS_CONVECTIVE, EQNS_DIVERGENCE, EQNS_SKEWSYMMETRIC
    use NavierStokes, only: visc, schmidt, prandtl
    use DNS_Arrays
    use BOUNDARY_BCS
    use Buffer

    implicit none
    private

    public :: TMarch_Initialize
    public :: TMarch_RungeKutta
    public :: TMarch_Courant

    real(wp), public :: dtime                   ! time step
    real(wp), public :: dte                     ! time step of each substep
    logical, public :: remove_divergence        ! Remove residual divergence every time step

    ! -------------------------------------------------------------------
    ! integer :: imode_rhs                        ! Type of implementation of the RHS of evolution equations
    ! integer, parameter :: EQNS_RHS_SPLIT = 18
    ! integer, parameter :: EQNS_RHS_COMBINED = 19
    ! integer, parameter :: EQNS_RHS_NONBLOCKING = 20

    ! type :: tmarch_dt
    !     sequence
    !     integer type
    !     integer nb_stages                           ! number of stages
    !     real(wp) kdt(5), kco(4), ktime(5)           ! explicit scheme coefficients
    !     real(wp) kex(3), kim(3)                     ! implicit scheme coefficients
    !     procedure(tmarch_interface), pointer, nopass :: tmarch_scheme
    ! end type tmarch_dt

    integer(wi) :: rkm_mode                     ! Type of Runge-Kutta scheme
    integer, parameter :: RKM_EXP3 = 3
    integer, parameter :: RKM_EXP4 = 4
    integer, parameter :: RKM_IMP3_DIFFUSION = 5
    integer, parameter :: RKM_IMP3_SOURCE = 6
    integer, parameter :: RKM_IMP3_DIFFSOURCE = 7

    integer(wi) :: rkm_endstep          ! number of substeps
    integer(wi) :: rkm_substep          ! substep counter

    real(wp) :: cfla, cfld, cflr        ! CFL numbers
    real(wp) etime                              ! time at each substep

    real(wp) kdt(5), kco(4), ktime(5)           ! explicit scheme coefficients
    real(wp) kex(3), kim(3)                     ! implicit scheme coefficients

    real(wp) schmidtfactor, dx2i
    integer(wi) i, j, k, kdsp, jdsp, idsp, is
    real(wp) dummy

    type :: ds_dt
        real(wp), allocatable :: one_ov_ds1(:)
        real(wp), allocatable :: one_ov_ds2(:)
    end type
    type(ds_dt) :: ds(3)

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_Initialize(inifile)
        use TLab_Memory, only: TLab_Allocate_Real
        use FDM, only: g

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, lstr
        character(len=128) eStr
        character(len=512) sRes
        integer ig

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Scheme=<RungeKuttaExplicit3/RungeKuttaExplicit4/RungeKuttaDiffusion3>')
        call TLab_Write_ASCII(bakfile, '#TimeStep=<value (used if CFL is negative)>')
        call TLab_Write_ASCII(bakfile, '#MaxCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#TimeDiffusiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#TimeReactiveCFL=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Scheme', 'dummy', sRes)
        if (trim(adjustl(sRes)) == 'rungekuttaexplicit3') then; rkm_mode = RKM_EXP3; lstr = '0.6'; 
        elseif (trim(adjustl(sRes)) == 'rungekuttaexplicit4') then; rkm_mode = RKM_EXP4; lstr = '1.2'; 
        elseif (trim(adjustl(sRes)) == 'rungekuttadiffusion3') then; rkm_mode = RKM_IMP3_DIFFUSION; lstr = '0.6'; 
            !  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rungekuttasource3'    ) THEN; rkm_mode = RKM_IMP3_SOURCE;
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong Scheme option.')
            call TLab_Stop(DNS_ERROR_RKORDER)
        end if

        ! Default cfla value set in lstr while reading Scheme
        call ScanFile_Real(bakfile, inifile, block, 'MaxCFL', trim(adjustl(lstr)), cfla)
        write (lstr, *) 0.25_wp*cfla ! Default value for diffusive CFL
        call ScanFile_Real(bakfile, inifile, block, 'TimeDiffusiveCFL', trim(adjustl(lstr)), cfld)
        write (lstr, *) 0.5_wp*cfla ! Default value for reactive CFL
        call ScanFile_Real(bakfile, inifile, block, 'TimeReactiveCFL', trim(adjustl(lstr)), cflr)
        call ScanFile_Real(bakfile, inifile, block, 'TimeStep', '0.05', dtime)

        ! ! -------------------------------------------------------------------
        ! ! Implicit RKM part
        ! ! -------------------------------------------------------------------
        ! if (rkm_mode == RKM_IMP3_DIFFUSION) then
        !     do is = 1, inb_scal
        !         if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN .or. &
        !             BcsScalJmax%type(is) == DNS_BCS_NEUMANN) then
        !             write (sRes, *) is; sRes = trim(adjustl(eStr))//'Scalar'//trim(adjustl(sRes))//'. Finite flux BC not implemented for SEMI-IMPLICITE DIFFUSION'
        !             call TLab_Write_ASCII(wfile, trim(adjustl(sRes)))
        !             write (sRes, *) is; sRes = trim(adjustl(eStr))//'Scalar'//trim(adjustl(sRes))//'. Setting fluxes at boundary to zero'
        !             call TLab_Write_ASCII(wfile, trim(adjustl(sRes)))
        !         end if
        !     end do

        ! end if

        ! -------------------------------------------------------------------
        call TLab_Write_ASCII(bakfile, '#RemoveDivergence=<none/remove>')

        call ScanFile_Char(bakfile, inifile, block, 'TermDivergence', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'no') then; remove_divergence = .false.
        else if (trim(adjustl(sRes)) == 'yes') then; remove_divergence = .true.
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong RemoveDivergence option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        ! RK coefficients
        select case (rkm_mode)
        case (RKM_EXP3)             ! Runge-Kutta explicit 3th order from Williamson 1980
            rkm_endstep = 3

            kdt(1:3) = [1.0_wp/3.0_wp, 15.0_wp/16.0_wp, 8.0_wp/15.0_wp]
            ktime(1:3) = [0.0_wp, 1.0_wp/3.0_wp, 3.0_wp/4.0_wp]
            kco(1:2) = [-5.0_wp/9.0_wp, -153.0_wp/128.0_wp]

        case (RKM_EXP4)             ! Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
            rkm_endstep = 5

            kdt(1) = 1432997174477.0_wp/9575080441755.0_wp  !C_1432997174477_R/C_9575080441755_R
            kdt(2) = 5161836677717.0_wp/13612068292357.0_wp !C_5161836677717_R/C_13612068292357_R
            kdt(3) = 1720146321549.0_wp/2090206949498.0_wp  !C_1720146321549_R/C_2090206949498_R
            kdt(4) = 3134564353537.0_wp/4481467310338.0_wp  !C_3134564353537_R/C_4481467310338_R
            kdt(5) = 2277821191437.0_wp/14882151754819.0_wp !C_2277821191437_R/C_14882151754819_R

            ktime(1) = 0.0_wp
            ktime(2) = kdt(1) !C_1432997174477_R/C_9575080441755_R
            ktime(3) = 2526269341429.0_wp/6820363962896.0_wp !C_2526269341429_R/C_6820363962896_R
            ktime(4) = 2006345519317.0_wp/3224310063776.0_wp !C_2006345519317_R/C_3224310063776_R
            ktime(5) = 2802321613138.0_wp/2924317926251.0_wp !C_2802321613138_R/C_2924317926251_R

            kco(1) = -567301805773.0_wp/1357537059087.0_wp     !C_567301805773_R/C_1357537059087_R
            kco(2) = -2404267990393.0_wp/2016746695238.0_wp    !C_2404267990393_R/C_2016746695238_R
            kco(3) = -3550918686646.0_wp/2091501179385.0_wp    !C_3550918686646_R/C_2091501179385_R
            kco(4) = -1275806237668.0_wp/842570457699.0_wp     !C_1275806237668_R/C_842570457699_R

        case (RKM_IMP3_DIFFUSION)   ! Runge-Kutta semi-implicit 3th order from Spalart, Moser & Rogers (1991)
            rkm_endstep = 3

            kdt(1:3) = [8.0_wp/15.0_wp, 5.0_wp/12.0_wp, 3.0_wp/4.0_wp]

            kim(1:3) = [111.0_wp/256.0_wp, 1.0_wp/2.0_wp, 2.0_wp/9.0_wp]
            kex(1:3) = [145.0_wp/256.0_wp, -9.0_wp/50.0_wp, 2.0_wp/9.0_wp]
            kco(1:3) = [0.0_wp, -17.0_wp/25.0_wp, -5.0_wp/9.0_wp]
            ! TO DO - calculate ktime from coefficients  ktime
            ktime(1:3) = [0.0_wp, 0.0_wp, 0.0_wp]

            ! Coefficients from Spalart, Moser, Rogers (1991)
            ! kim = beta/gamma
            ! kex = alpha/gamma
            ! kco = zeta/gamma
            !
            ! alpha = [ 29./96.,   -3./40,    1./6. ]
            ! beta  = [ 37./160.,   5./24.,   1./6. ]
            ! gamma = [  8./15.,    5./12.,   3./4. ]
            ! zeta  = [  0.,      -17./60.,  -5./12.]

        end select

        ! ###################################################################
        ! Memory management
        call TLab_Allocate_Real(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
        call TLab_Allocate_Real(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

        p_hq(1:imax, 1:jmax, 1:kmax, 1:inb_flow) => hq(1:imax*jmax*kmax*inb_flow, 1)
        p_hs(1:imax, 1:jmax, 1:kmax, 1:inb_scal) => hs(1:imax*jmax*kmax*inb_scal, 1)

        ! ###################################################################
        ! maximum diffusivities for TMarch_Courant
        schmidtfactor = 1.0_wp
        dummy = 1.0_wp/prandtl
        schmidtfactor = max(schmidtfactor, dummy)
        dummy = 1.0_wp/minval(schmidt(1:inb_scal))
        schmidtfactor = max(schmidtfactor, dummy)
        schmidtfactor = schmidtfactor*visc

        ! ###################################################################
        do ig = 1, 3
            allocate (ds(ig)%one_ov_ds1(g(ig)%size))
            ds(ig)%one_ov_ds1(:) = 1.0_wp/g(ig)%jac(:, 1)

            allocate (ds(ig)%one_ov_ds2(g(ig)%size))
            ds(ig)%one_ov_ds2(:) = ds(ig)%one_ov_ds1(:)*ds(ig)%one_ov_ds1(:)

        end do

        ! Maximum of (1/dx^2 + 1/dy^2 + 1/dz^2) for TMarch_Courant
#ifdef USE_MPI
        idsp = ims_offset_i
        jdsp = ims_offset_j
        kdsp = ims_offset_k
#else
        idsp = 0
        jdsp = 0
        kdsp = 0
#endif

        dx2i = 0.0_wp
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    dummy = 0.0_wp
                    if (x%size > 1) dummy = dummy + ds(1)%one_ov_ds2(i + idsp)
                    if (y%size > 1) dummy = dummy + ds(2)%one_ov_ds2(j + jdsp)
                    if (z%size > 1) dummy = dummy + ds(3)%one_ov_ds2(k + kdsp)
                    dx2i = max(dx2i, dummy)
                end do
            end do
        end do

        return
    end subroutine TMarch_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_RungeKutta()
        use TLab_Arrays
        use DNS_LOCAL
        use DNS_Control
        use DNS_Arrays

        ! -------------------------------------------------------------------
        real(wp) alpha

#ifdef USE_PROFILE
        integer(wi) t_srt, t_end, t_dif, idummy, PROC_CYCLES, MAX_CYCLES
        character*256 time_string
#endif

        !########################################################################
        ! -------------------------------------------------------------------
        ! Initialize arrays to zero for the explcit low-storage algorithm
        ! -------------------------------------------------------------------
        if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
            if (flow_on) hq = 0.0_wp
            if (scal_on) hs = 0.0_wp
            ! if (part%type /= PART_TYPE_NONE) l_hq = 0.0_wp
        end if
        !########################################################################
        ! Loop over the sub-stages
        !########################################################################
        do rkm_substep = 1, rkm_endstep

            ! -------------------------------------------------------------------
            ! Update transported (or prognostic) variables q and s
            ! -------------------------------------------------------------------
            dte = dtime*kdt(rkm_substep)
            etime = rtime + dtime*ktime(rkm_substep)

#ifdef USE_PROFILE
            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
#endif
            ! I could define procedure pointers to handle this...
            select case (nse_eqns)
            case (DNS_EQNS_BOUSSINESQ)
                if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                    call TMarch_Substep_Boussinesq_Explicit()
                else if (rkm_mode == RKM_IMP3_DIFFUSION) then
                    call TMarch_Substep_Boussinesq_Implicit()
                end if

            case (DNS_EQNS_ANELASTIC)
                if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                    call TMarch_Substep_Anelastic_Explicit()
                    ! else if (rkm_mode == RKM_IMP3_DIFFUSION) then
                    !     call TMarch_Substep_Boussinesq_Implicit()
                end if

            case (DNS_EQNS_COMPRESSIBLE)

            end select

            call DNS_BOUNDS_LIMIT()
            if (int(logs_data(1)) /= 0) return ! Error detected

            ! -------------------------------------------------------------------
            ! Update RHS hq and hs in the explicit low-storage algorithm
            ! -------------------------------------------------------------------
            if ((rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) .and. &
                rkm_substep < rkm_endstep) then

                alpha = kco(rkm_substep)

                if (flow_on) then
                    do is = 1, inb_flow
#ifdef USE_BLAS
                        call DSCAL(isize_field, alpha, hq(:, is), 1)
#else
                        hq(:, is) = alpha*hq(:, is)
#endif
                    end do
                end if

                if (scal_on) then
                    do is = 1, inb_scal
#ifdef USE_BLAS
                        call DSCAL(isize_field, alpha, hs(:, is), 1)
#else
                        hs(:, is) = alpha*hs(:, is)
#endif
                    end do
                end if

            end if

            ! -------------------------------------------------------------------
            ! Profiling data
            ! -------------------------------------------------------------------
#ifdef USE_PROFILE
            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)
            idummy = t_end - t_srt

#ifdef USE_MPI
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
            if (ims_pro == 0) then
                write (time_string, 999) ims_npro, ims_npro_i, ims_npro_k, rkm_substep, t_dif/1.0_wp/PROC_CYCLES/ims_npro
999             format(I5.5, ' (ims_npro_i X ims_npro_k:', I4.4, 'x', I4.4, 1x, ') RK-Substep', I1, ':', E13.5, 's')
                call TLab_Write_ASCII(lfile, time_string)
            end if
#else
            t_dif = idummy
            write (time_string, 999) rkm_substep, t_dif/1.0_wp/PROC_CYCLES/ims_npro
999         format('RK-Substep', I1, ':', E13.5, 's')
            call TLab_Write_ASCII(lfile, time_string)
#endif

#endif
        end do

        return
    end subroutine TMarch_RungeKutta

    !########################################################################
    !#
    !# Determine the variable time step is a positive cfl is given
    !# If negative cfl is prescribed, then constant dtime and this routine
    !# calculates CFL and diffustion numbers for log files
    !#
    !# The reacting case should be reviewed in the new formulation in terms
    !# of energy.
    !#
    !# The diffusion number is fixed in terms of the CFL, which is the input.
    !# This depends on the scheme used. From Lele (1992), page 32, we have that, if
    !# the sixth order tridiagonal scheme is used, then the maximum CFL number
    !# for a 4RK is 2.9/1.989, about 1.43. For the (5)4RK from CarpenterKennedy1994
    !# used here we have 3.36/1.989, about 1.69.
    !# This holds for periodic case. A safety margin leads to the common value of 1.2.
    !#
    !# If second order finite different operator is used, then the maximum
    !# diffusion number is 2.9/6.857, about 0.42.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/6.857 = 0.68
    !# If the extension by Lamballais et al is used, then the maximum
    !# diffusion number is 2.9/pi^2, about 0.29.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/pi^2 = 0.47.
    !#
    !# If twice the first order finite difference operator is used, then the
    !# maximum diffusion number is 2.9/1.989^2, about 0.73.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
    !#
    !# In incompressible mode the arrays rho, p and vis are not used
    !#
    !########################################################################
    subroutine TMarch_Courant()
        use DNS_Control, only: logs_data, logs_dtime
        ! use Thermodynamics, only: gamma0, itransport, EQNS_TRANS_POWERLAW
        use TLab_Pointers_3D, only: u, v, w, p_wrk3d !, p, rho, vis

        ! -------------------------------------------------------------------
        integer(wi) ipmax, j_glo
        real(wp) dt_loc
        real(wp) pmax(3), dtc, dtd, dtr
#ifdef USE_MPI
        real(wp) pmax_aux(3)
#endif

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i
        jdsp = ims_offset_j
#else
        idsp = 0
        jdsp = 0
#endif

        dtc = big_wp   ! So that the minimum non-zero determines dt at the end
        dtd = big_wp
        dtr = big_wp

        ipmax = 0       ! Initialize counter of time constraints

        ! ###################################################################
        ! CFL number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of u/dx + v/dy + w/dz
        ! -------------------------------------------------------------------
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            if (y%size > 1) then
                do k = 1, kmax
                    do j = 1, jmax
                        j_glo = j + jdsp
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + idsp) &
                                               + abs(v(i, j, k))*ds(2)%one_ov_ds1(j_glo) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            else    ! do I need this?
                do k = 1, kmax
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + idsp) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            end if

        end select

        pmax(1) = maxval(p_wrk3d)

        ! ###################################################################
        ! Diffusion number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of \mu*(1/dx^2 + 1/dy^2 + 1/dz^2)
        ! -------------------------------------------------------------------
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            pmax(2) = schmidtfactor*dx2i

        end select

        ! ###################################################################
        ! Final operations
        ! ###################################################################
#ifdef USE_MPI
        call MPI_ALLREDUCE(pmax, pmax_aux, ipmax, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        pmax(1:ipmax) = pmax_aux(1:ipmax)
#endif

        if (pmax(1) > 0.0_wp) dtc = cfla/pmax(1) ! Set time step for the given CFL number
        if (pmax(2) > 0.0_wp) dtd = cfld/pmax(2) ! Set time step for the given diffusion number

        ! -------------------------------------------------------------------
        if (cfla > 0.0_wp) then
            if (rkm_mode == RKM_EXP3 .or. rkm_mode == RKM_EXP4) then
                dt_loc = min(dtc, dtd)
            else
                dt_loc = dtc
            end if
            dt_loc = min(dt_loc, dtr)

            dtime = dt_loc

        end if

        ! Real CFL and diffusion numbers being used, for the logfile
        logs_dtime = dtime
        logs_data(2) = dtime*pmax(1)
        logs_data(3) = dtime*pmax(2)

        return

    end subroutine TMarch_Courant

    !########################################################################
    !########################################################################
    subroutine TMarch_Substep_Boussinesq_Explicit()
        use TLab_Arrays, only: q, s, txc
        use DNS_Arrays, only: hq, hs
        use TLab_Sources

        ! #######################################################################
        ! Accumulate RHS terms
        call TLab_Sources_Flow(q, s, hq, txc(:, 1))
        call TLab_Sources_Scal(s, hs, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))

        if (bufferType == BUFFER_TYPE_NUDGE) call Buffer_Nudge()

        call NSE_Boussinesq()

        ! #######################################################################
        ! Perform the time stepping
        do is = 1, inb_flow
            q(:, is) = q(:, is) + dte*hq(:, is)
        end do

        do is = 1, inb_scal
            s(:, is) = s(:, is) + dte*hs(:, is)
        end do

        return
    end subroutine TMarch_Substep_Boussinesq_Explicit

    !########################################################################
    !########################################################################
    subroutine TMarch_Substep_Anelastic_Explicit()
        use TLab_Arrays, only: q, s, txc
        use DNS_Arrays, only: hq, hs
        use TLab_Sources

        ! #######################################################################
        ! Accumulate RHS terms
        call TLab_Sources_Flow(q, s, hq, txc(:, 1))
        call TLab_Sources_Scal(s, hs, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))

        if (bufferType == BUFFER_TYPE_NUDGE) call Buffer_Nudge()

        call NSE_Anelastic()

        ! #######################################################################
        ! Perform the time stepping
        do is = 1, inb_flow
            q(:, is) = q(:, is) + dte*hq(:, is)
        end do

        do is = 1, inb_scal
            s(:, is) = s(:, is) + dte*hs(:, is)
        end do

        call TLab_Diagnostic(imax, jmax, kmax, s)

        return
    end subroutine TMarch_Substep_Anelastic_Explicit

    !########################################################################
    !########################################################################
    subroutine TMarch_Substep_Boussinesq_Implicit()

        ! call RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2(kex(rkm_substep), kim(rkm_substep), kco(rkm_substep))

        ! ! ! pressure-correction algorithm; to be checked
        ! ! CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3(&
        ! !      kex,kim,kco,  &
        ! !      q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
        ! !      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), txc(1,8))

        return
    end subroutine TMarch_Substep_Boussinesq_Implicit

    ! !########################################################################
    ! !########################################################################
    !     subroutine TMarch_SUBSTEP_COMPRESSIBLE()
    !         use TLab_Arrays
    !         use TLab_Pointers
    !         use THERMO_CALORIC, only: THERMO_GAMMA
    !         use Thermodynamics, only: CRATIO_INV
    !         use DNS_Arrays
    !         use BOUNDARY_BUFFER
    !         use BOUNDARY_BCS, only: BcsDrift
    !         use BOUNDARY_BCS_COMPRESSIBLE

    !         ! -------------------------------------------------------------------
    !         real(wp) rho_ratio, dt_rho_ratio, prefactor, M2_max
    !         integer(wi) inb_scal_loc

    !         ! ###################################################################
    !         ! Evaluate standard RHS of equations
    !         ! global formulation
    !         ! ###################################################################
    !         if (nse_eqns == DNS_EQNS_COMPRESSIBLE .and. &
    !             nse_advection == EQNS_SKEWSYMMETRIC .and. &
    !             nse_viscous == EQNS_EXPLICIT .and. &
    !             nse_diffusion == EQNS_EXPLICIT) then
    !             call RHS_FLOW_GLOBAL_2()

    !             do is = 1, inb_scal
    !                 call RHS_SCAL_GLOBAL_2(is)
    !             end do

    !         else
    !             ! ###################################################################
    !             ! Evaluate standard RHS of equations
    !             ! split formulation
    !             ! ###################################################################
    !             if (nse_advection == EQNS_DIVERGENCE) then
    !                 call RHS_FLOW_EULER_DIVERGENCE()
    !                 do is = 1, inb_scal
    !                     call RHS_SCAL_EULER_DIVERGENCE()
    !                 end do

    !             else if (nse_advection == EQNS_SKEWSYMMETRIC) then
    !                 call RHS_FLOW_EULER_SKEWSYMMETRIC()
    !                 do is = 1, inb_scal
    !                     call RHS_SCAL_EULER_SKEWSYMMETRIC(is)
    !                 end do
    !             end if

    !             if (nse_viscous == EQNS_DIVERGENCE) then
    !                 call RHS_FLOW_VISCOUS_DIVERGENCE()

    !             else if (nse_viscous == EQNS_EXPLICIT) then
    !                 call RHS_FLOW_VISCOUS_EXPLICIT()

    !             end if

    !             if (nse_diffusion == EQNS_DIVERGENCE) then
    !                 ! diffusion transport of enthalpy is accumulated in txc5:7 and used in RHS_FLOW_CONDUCTION
    !                 txc(:, 5:7) = 0.0_wp
    !                 do is = 1, inb_scal
    !                     ! Check this routine for the airwater case
    !                     call RHS_SCAL_DIFFUSION_DIVERGENCE(is)
    !                 end do
    !                 call RHS_FLOW_CONDUCTION_DIVERGENCE()

    !             else if (nse_diffusion == EQNS_EXPLICIT) then
    !                 do is = 1, inb_scal
    !                     call RHS_SCAL_DIFFUSION_EXPLICIT(is)
    !                 end do
    !                 call RHS_FLOW_CONDUCTION_EXPLICIT()

    !             end if

    !         end if

    !         ! ###################################################################
    !         ! Impose boundary conditions
    !         ! Temperature array T is used as auxiliary array because it is no
    !         ! longer used until the fields are updated
    !         ! ###################################################################
    ! #define GAMMA_LOC(i) txc(i,6)
    ! #define AUX_LOC(i)   T(i)

    !         call THERMO_GAMMA(imax*jmax*kmax, s, T, GAMMA_LOC(:))

    !         ! Maximum Mach for Poinsot & Lele reference pressure BC
    !         if (BcsDrift) then
    !             M2_max = 0.0_wp
    !             do i = 1, isize_field
    !                 dummy = (u(i)*u(i) + v(i)*v(i) + w(i)*w(i))*rho(i)/(GAMMA_LOC(i)*p(i))
    !                 M2_max = max(M2_max, dummy)
    !             end do
    ! #ifdef USE_MPI
    !             call MPI_ALLREDUCE(M2_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
    !             M2_max = dummy
    ! #endif
    !         end if

    !         if (.not. y%periodic) then
    !             call BOUNDARY_BCS_Y(isize_field, M2_max, rho, u, v, w, p, GAMMA_LOC(1), s, &
    !                                 hq(1, 5), hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), hs, &
    !                                 txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), AUX_LOC(:))
    !         end if

    !         if (.not. x%periodic) then
    !             call BOUNDARY_BCS_X(isize_field, M2_max, etime, rho, u, v, w, p, GAMMA_LOC(1), s, &
    !                                 hq(1, 5), hq(1, 1), hq(1, 2), hq(1, 3), hq(1, 4), hs, txc, AUX_LOC(:))
    !         end if

    ! #undef GAMMA_LOC
    ! #undef AUX_LOC

    !         ! ###################################################################
    !         ! Impose buffer zone as relaxation terms
    !         ! ###################################################################
    !         if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
    !             call BOUNDARY_BUFFER_RELAX_FLOW()
    !             call BOUNDARY_BUFFER_RELAX_SCAL()
    !         end if

    !         ! ###################################################################
    !         ! Perform the time stepping
    !         ! ###################################################################
    !         rho_ratio = 1.0_wp
    !         prefactor = CRATIO_INV*0.5_wp

    !         if (flow_on) then
    !             if (scal_on) then; inb_scal_loc = inb_scal
    !             else; inb_scal_loc = 0
    !             end if

    !             ! -------------------------------------------------------------------
    !             ! Total energy formulation
    !             ! -------------------------------------------------------------------
    !             if (nse_eqns == DNS_EQNS_TOTAL) then
    ! !$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
    ! !$omp do
    !                 do i = 1, isize_field
    !                     rho_ratio = rho(i)
    !                     rho(i) = rho(i) + dte*hq(i, 5)
    !                     rho_ratio = rho_ratio/rho(i)
    !                     dt_rho_ratio = dte/rho(i)

    !                     e(i) = rho_ratio*(e(i) + prefactor*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))) &
    !                            + dt_rho_ratio*hq(i, 4)

    !                     u(i) = rho_ratio*u(i) + dt_rho_ratio*hq(i, 1)
    !                     v(i) = rho_ratio*v(i) + dt_rho_ratio*hq(i, 2)
    !                     w(i) = rho_ratio*w(i) + dt_rho_ratio*hq(i, 3)

    !                     e(i) = e(i) - prefactor*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))

    !                     do is = 1, inb_scal_loc
    !                         s(i, is) = rho_ratio*s(i, is) + dt_rho_ratio*hs(i, is)
    !                     end do
    !                 end do
    ! !$omp end do
    ! !$omp end parallel

    !                 ! -------------------------------------------------------------------
    !                 ! Internal energy formulation
    !                 ! -------------------------------------------------------------------
    !             else if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
    ! !$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
    ! !$omp do
    !                 do i = 1, isize_field
    !                     rho_ratio = rho(i)
    !                     rho(i) = rho(i) + dte*hq(i, 5)
    !                     rho_ratio = rho_ratio/rho(i)
    !                     dt_rho_ratio = dte/rho(i)

    !                     e(i) = rho_ratio*e(i) + dt_rho_ratio*hq(i, 4)
    !                     u(i) = rho_ratio*u(i) + dt_rho_ratio*hq(i, 1)
    !                     v(i) = rho_ratio*v(i) + dt_rho_ratio*hq(i, 2)
    !                     w(i) = rho_ratio*w(i) + dt_rho_ratio*hq(i, 3)

    !                     do is = 1, inb_scal_loc
    !                         s(i, is) = rho_ratio*s(i, is) + dt_rho_ratio*hs(i, is)
    !                     end do
    !                 end do
    ! !$omp end do
    ! !$omp end parallel

    !             end if

    !         else
    !             if (scal_on) then
    !                 do is = 1, inb_scal
    ! !$omp parallel default( shared ) private( i, dt_rho_ratio )
    ! !$omp do
    !                     do i = 1, isize_field
    !                         dt_rho_ratio = dte/rho(i)
    !                         s(i, is) = rho_ratio*s(i, is) + dt_rho_ratio*hs(i, is)
    !                     end do
    ! !$omp end do
    ! !$omp end parallel
    !                 end do
    !             end if
    !         end if

    !         ! ###################################################################
    !         ! Impose buffer zone as filter
    !         ! ###################################################################
    !         if (BuffType == DNS_BUFFER_FILTER .or. BuffType == DNS_BUFFER_BOTH) then
    !             call BOUNDARY_BUFFER_FILTER(rho, u, v, w, e, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
    !         end if

    !         return
    !     end subroutine TMarch_SUBSTEP_COMPRESSIBLE

    !     !########################################################################
    !     !########################################################################
    !     subroutine TMarch_SUBSTEP_PARTICLE()
    !         use DNS_Arrays, only: l_hq
    !         use PARTICLE_VARS
    !         use PARTICLE_ARRAYS

    !         ! -------------------------------------------------------------------
    !         integer(wi) is

    ! #ifdef USE_MPI
    !         integer(wi) nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north
    ! #endif
    !         real(wp) x_right, z_right
    !         real(wp) y_right

    !         !#####################################################################
    !         call RHS_PART_1()

    !         !#######################################################################
    !         ! Update particle properties
    !         !#######################################################################
    !         do is = 1, inb_part
    !             l_q(1:l_g%np, is) = l_q(1:l_g%np, is) + dte*l_hq(1:l_g%np, is)
    !         end do

    !         !#####################################################################
    !         ! Boundary control to see if particles leave processor
    !         !#####################################################################
    ! #ifdef USE_MPI

    !         ! -------------------------------------------------------------------
    !         ! Particle sorting for Send/Recv X-Direction
    !         ! -------------------------------------------------------------------
    !         if (ims_npro_i > 1) then
    !             call PARTICLE_MPI_SORT(1, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

    !             if (ims_pro_i == 0) then !Take care of periodic boundary conditions west
    !                 if (nzone_west /= 0) then
    !                     l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) = &
    !                         l_q(nzone_grid + 1:nzone_grid + nzone_west, 1) + x%scale
    !                 end if
    !             end if

    !             if (ims_pro_i == (ims_npro_i - 1)) then !Take care of periodic boundary conditions east
    !                 if (nzone_east /= 0) then
    !                     l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) = &
    !                         l_q(nzone_grid + nzone_west + 1:nzone_grid + nzone_west + nzone_east, 1) - x%scale
    !                 end if
    !             end if

    !             call PARTICLE_MPI_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, l_q, l_hq, l_g%tags, l_g%np)

    !         else
    !             x_right = x%nodes(1) + x%scale

    !             do i = 1, l_g%np
    !                 if (l_q(i, 1) > x_right) then
    !                     l_q(i, 1) = l_q(i, 1) - x%scale

    !                 elseif (l_q(i, 1) < x%nodes(1)) then
    !                     l_q(i, 1) = l_q(i, 1) + x%scale

    !                 end if

    !             end do

    !         end if

    !         ! -------------------------------------------------------------------
    !         ! Particle sorting for Send/Recv Z-Direction
    !         ! -------------------------------------------------------------------
    !         if (ims_npro_k > 1) then
    !             call PARTICLE_MPI_SORT(3, l_g, l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

    !             if (ims_pro_k == 0) then !Take care of periodic boundary conditions south
    !                 if (nzone_south /= 0) then
    !                     l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) = &
    !                         l_q(nzone_grid + 1:nzone_grid + nzone_south, 3) + z%scale
    !                 end if
    !             end if

    !             if (ims_pro_k == (ims_npro_k - 1)) then !Take care of periodic boundary conditions north
    !                 if (nzone_north /= 0) then
    !                     l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) = &
    !                         l_q(nzone_grid + nzone_south + 1:nzone_grid + nzone_south + nzone_north, 3) - z%scale
    !                 end if
    !             end if

    !             call PARTICLE_MPI_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, l_q, l_hq, l_g%tags, l_g%np)

    !         else
    !             call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    !             z_right = z%nodes(1) + z%scale

    !             do i = 1, l_g%np
    !                 if (l_q(i, 3) > z_right) then
    !                     l_q(i, 3) = l_q(i, 3) - z%scale

    !                 elseif (l_q(i, 3) < z%nodes(1)) then
    !                     l_q(i, 3) = l_q(i, 3) + z%scale

    !                 end if

    !             end do

    !         end if

    ! #else
    !         !#######################################################################
    !         ! Serial; would it be faster to use MOD?
    !         x_right = x%nodes(1) + x%scale
    !         z_right = z%nodes(1) + z%scale

    !         do i = 1, l_g%np
    !             if (l_q(i, 1) > x_right) then
    !                 l_q(i, 1) = l_q(i, 1) - x%scale

    !             elseif (l_q(i, 1) < x%nodes(1)) then
    !                 l_q(i, 1) = l_q(i, 1) + x%scale

    !             end if

    !             if (l_q(i, 3) > z_right) then
    !                 l_q(i, 3) = l_q(i, 3) - z%scale

    !             elseif (l_q(i, 3) < z%nodes(1)) then
    !                 l_q(i, 3) = l_q(i, 3) + z%scale

    !             end if

    !         end do

    ! #endif

    !         y_right = y%nodes(1) + y%scale
    !         select case (part_bcs)
    !         case (PART_BCS_SPECULAR)
    !             do i = 1, l_g%np
    !                 if (l_q(i, 2) > y_right) then
    !                     l_q(i, 2) = 2*y_right - l_q(i, 2)
    !                     l_q(i, 5) = -l_q(i, 5)
    !                 elseif (l_q(i, 2) < y%nodes(1)) then
    !                     l_q(i, 2) = 2*y%nodes(1) - l_q(i, 2)
    !                     l_q(i, 5) = -l_q(i, 5)
    !                 end if
    !             end do

    !         case (PART_BCS_STICK)
    !             do i = 1, l_g%np
    !                 if (l_q(i, 2) > y_right) then
    !                     l_q(i, 2) = y_right
    !                 elseif (l_q(i, 2) < y%nodes(1)) then
    !                     l_q(i, 2) = y%nodes(1)
    !                 end if
    !             end do

    !         end select

    !         !#######################################################################
    !         ! Recalculating closest node below in Y direction
    !         !#######################################################################
    !         call LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, y%size, y%nodes)

    !         return
    !     end subroutine TMarch_SUBSTEP_PARTICLE

end module TimeMarching
