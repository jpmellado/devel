!########################################################################
!#
!# Evolution equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and w transposes are calculated first for the
!# Ox and Oz momentum equations, stored in tmp5 and tmp6 and then used as needed.
!# This saves 2 MPI transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
subroutine RHS_GLOBAL_INCOMPRESSIBLE_1()
#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLab_OpenMP
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_flow, inb_scal
    use TLab_Arrays
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
    use FDM, only: g
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    ! use Thermo_Anelastic
    use DNS_Arrays
    use DNS_LOCAL, only: remove_divergence
    use TimeMarching, only: dte
    ! use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    ! use IBM_VARS, only: imode_ibm, imode_ibm_scal, ibm_burgers
    use OPR_Partial
    use OPR_Burgers
    use OPR_Elliptic
    ! use OPR_FILTERS

    implicit none

    ! -----------------------------------------------------------------------
    integer(wi) iq, is, ij
    integer ibc
    real(wp) dummy

    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning

#ifdef USE_ESSL
    integer ilen
#endif

    ! #######################################################################
#ifdef USE_ESSL
    ilen = isize_field
#endif

    ! #######################################################################
    ! Preliminaries for Scalar BC
    ! (flow BCs initialized below as they are used for pressure in between)
    ! #######################################################################
    ! Default is zero
    ! BcsScalKmin%ref(:, :, :) = 0.0_wp
    ! BcsScalKmax%ref(:, :, :) = 0.0_wp

    ! ! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    ! if (any(BcsScalJmin%SfcType(1:inb_scal) == DNS_SFC_LINEAR) .or. &
    !     any(BcsScalJmax%SfcType(1:inb_scal) == DNS_SFC_LINEAR)) then
    !     do is = 1, inb_scal
    !         p_bcs(1:imax, 1:jmax, 1:kmax) => hs(1:imax*jmax*kmax, is)
    !         if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) BcsScalKmin%ref(:, :, is) = p_bcs(:, 1, :)
    !         if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) BcsScalKmax%ref(:, :, is) = p_bcs(:, jmax, :)
    !     end do
    ! end if

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    ! Diagonal terms and transposed velocity arrays
    call OPR_Burgers_X(OPR_B_SELF, 0, imax, jmax, kmax, u, u, tmp1, tmp4) ! store u transposed in tmp4
    call OPR_Burgers_Y(OPR_B_SELF, 0, imax, jmax, kmax, v, v, tmp2, tmp5) ! store v transposed in tmp5
    call OPR_Burgers_Z(OPR_B_SELF, 0, imax, jmax, kmax, w, w, tmp3, tmp6) ! store w transposed in tmp6

    ! Ox momentum equation
    call OPR_Burgers_Y(OPR_B_U_IN, 0, imax, jmax, kmax, u, v, tmp7, tmp9, tmp5) ! tmp5 contains v transposed
    call OPR_Burgers_Z(OPR_B_U_IN, 0, imax, jmax, kmax, u, w, tmp8, tmp9, tmp6) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 1) = hq(ij, 1) + tmp1(ij) + tmp7(ij) + tmp8(ij)
    end do
!$omp end parallel

    ! Oy momentum equation
    call OPR_Burgers_X(OPR_B_U_IN, 0, imax, jmax, kmax, v, u, tmp7, tmp9, tmp4) ! tmp4 contains u transposed
    call OPR_Burgers_Z(OPR_B_U_IN, 0, imax, jmax, kmax, v, w, tmp8, tmp9, tmp6) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 2) = hq(ij, 2) + tmp2(ij) + tmp7(ij) + tmp8(ij)
    end do
!$omp end parallel

    ! Oz momentum equation
    call OPR_Burgers_X(OPR_B_U_IN, 0, imax, jmax, kmax, w, u, tmp7, tmp9, tmp4) ! tmp4 contains u transposed
    call OPR_Burgers_Y(OPR_B_U_IN, 0, imax, jmax, kmax, w, v, tmp8, tmp9, tmp5) ! tmp5 contains v transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 3) = hq(ij, 3) + tmp3(ij) + tmp7(ij) + tmp8(ij)
    end do
!$omp end parallel

    ! Scalar equations
    do is = 1, inb_scal
        call OPR_Burgers_X(OPR_B_U_IN, is, imax, jmax, kmax, s(1, is), u, tmp1, tmp9, tmp4) ! tmp4 contains u transposed
        call OPR_Burgers_Y(OPR_B_U_IN, is, imax, jmax, kmax, s(1, is), v, tmp2, tmp9, tmp5) ! tmp5 contains v transposed
        call OPR_Burgers_Z(OPR_B_U_IN, is, imax, jmax, kmax, s(1, is), w, tmp3, tmp9, tmp6) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        call TLab_OMP_PARTITION(isize_field, srt, end, siz)
        do ij = srt, end
            hs(ij, is) = hs(ij, is) + tmp1(ij) + tmp2(ij) + tmp3(ij)
        end do
!$omp end parallel

    end do

    ! #######################################################################
    ! Impose buffer zone as relaxation terms
    ! #######################################################################
    ! if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
    !     call BOUNDARY_BUFFER_RELAX_FLOW()
    ! end if

    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    if (remove_divergence) then ! remove residual divergence

#ifdef USE_ESSL
!$omp parallel default( shared )&
!$omp private( ilen, dummy, srt,end,siz )
#else
!$omp parallel default( shared )&
!$omp private( ij,   dummy, srt,end,siz )
#endif

        call TLab_OMP_PARTITION(isize_field, srt, end, siz)
        dummy = 1.0_wp/dte

#ifdef USE_ESSL
        ilen = siz
        call DZAXPY(ilen, dummy, v(srt), 1, hq(srt, 2), 1, tmp2(srt), 1)
        call DZAXPY(ilen, dummy, u(srt), 1, hq(srt, 1), 1, tmp3(srt), 1)
        call DZAXPY(ilen, dummy, w(srt), 1, hq(srt, 3), 1, tmp4(srt), 1)

#else
        do ij = srt, end
            tmp2(ij) = hq(ij, 2) + v(ij)*dummy
            tmp3(ij) = hq(ij, 1) + u(ij)*dummy
            tmp4(ij) = hq(ij, 3) + w(ij)*dummy
        end do

#endif
!$omp end parallel

        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            ! call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp2)
            ! call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp3)
            ! call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp4)
        end if
        ! if (stagger_on) then ! staggering on horizontal pressure nodes
        !     !  Oy derivative
        !     call OPR_Partial_X(OPR_P0_INT_VP, imax, jmax, kmax, g(1), tmp2, tmp5)
        !     call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), tmp5, tmp2)
        !     call OPR_Partial_Z(OPR_P0_INT_VP, imax, jmax, kmax, g(3), tmp2, tmp1)
        !     !  Ox derivative
        !     call OPR_Partial_X(OPR_P1_INT_VP, imax, jmax, kmax, g(1), tmp3, tmp5)
        !     call OPR_Partial_Z(OPR_P0_INT_VP, imax, jmax, kmax, g(3), tmp5, tmp2)
        !     !  Oz derivative
        !     call OPR_Partial_X(OPR_P0_INT_VP, imax, jmax, kmax, g(1), tmp4, tmp5)
        !     call OPR_Partial_Z(OPR_P1_INT_VP, imax, jmax, kmax, g(3), tmp5, tmp3)
        ! else
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), tmp2, tmp1)
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), tmp3, tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), tmp4, tmp3)
        ! end if

    else
        ! if (nse_eqns == DNS_EQNS_ANELASTIC) then
        !     call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 2), tmp2)
        !     call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 1), tmp3)
        !     call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 3), tmp4)
        !     call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), tmp2, tmp1)
        !     call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), tmp3, tmp2)
        !     call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), tmp4, tmp3)
        ! else
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), hq(:, 2), tmp1)
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), hq(:, 1), tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), hq(:, 3), tmp3)
        ! end if

    end if

    ! -----------------------------------------------------------------------
!$omp parallel default( shared ) private( ij,srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
    end do
!$omp end parallel

    ! -----------------------------------------------------------------------
    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)

    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        ! BcsFlowKmin%ref(:, :, 3) = p_hq(:, :, 1, 3)*rbackground(1)
        ! BcsFlowKmax%ref(:, :, 3) = p_hq(:, :, kmax, 3)*rbackground(kmax)
    else
        BcsFlowKmin%ref(:, :, 3) = p_hq(:, :, 1, 3)
        BcsFlowKmax%ref(:, :, 3) = p_hq(:, :, kmax, 3)
    end if

    ! pressure in tmp1, Oz derivative in tmp3
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp4, BcsFlowKmin%ref(:, :, 3), BcsFlowKmax%ref(:, :, 3))
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), tmp1, tmp3)

    ! ! filter pressure p and its vertical gradient dpdy
    ! if (any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
    !     call OPR_FILTER(imax, jmax, kmax, PressureFilter, tmp1, txc(1:isize_field,4:6))
    !     call OPR_FILTER(imax, jmax, kmax, PressureFilter, tmp3, txc(1:isize_field,4:6))
    ! end if

    ! if (stagger_on) then
    !     !  vertical pressure derivative   dpdy - back on horizontal velocity nodes
    !     call OPR_Partial_Z(OPR_P0_INT_PV, imax, jmax, kmax, g(3), tmp3, tmp5)
    !     call OPR_Partial_X(OPR_P0_INT_PV, imax, jmax, kmax, g(1), tmp5, tmp3)
    !     !  horizontal pressure derivative dpdz - back on horizontal velocity nodes
    !     call OPR_Partial_Z(OPR_P1_INT_PV, imax, jmax, kmax, g(3), tmp1, tmp5)
    !     call OPR_Partial_X(OPR_P0_INT_PV, imax, jmax, kmax, g(1), tmp5, tmp4)
    !     !  horizontal pressure derivative dpdx - back on horizontal velocity nodes
    !     call OPR_Partial_Z(OPR_P0_INT_PV, imax, jmax, kmax, g(3), tmp1, tmp5)
    !     call OPR_Partial_X(OPR_P1_INT_PV, imax, jmax, kmax, g(1), tmp5, tmp2)
    ! else
    !  horizontal pressure derivatives
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), tmp1, tmp2)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), tmp1, tmp4)
    ! end if

    ! -----------------------------------------------------------------------
    ! Add pressure gradient
    ! -----------------------------------------------------------------------
    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        ! call Thermo_Anelastic_WEIGHT_SUBTRACT(imax, jmax, kmax, ribackground, tmp2, hq(:, 1))
        ! call Thermo_Anelastic_WEIGHT_SUBTRACT(imax, jmax, kmax, ribackground, tmp3, hq(:, 2))
        ! call Thermo_Anelastic_WEIGHT_SUBTRACT(imax, jmax, kmax, ribackground, tmp4, hq(:, 3))

    else
#ifdef USE_ESSL
!$omp parallel default( shared ) &
!$omp private( ilen, srt,end,siz,dummy )
#else
!$omp parallel default( shared ) &
!$omp private( ij,   srt,end,siz,dummy )
#endif
        call TLab_OMP_PARTITION(isize_field, srt, end, siz)

#ifdef USE_ESSL
        ilen = siz
        dummy = -1.0_wp
        call DAXPY(ilen, dummy, tmp2(srt), 1, hq(srt, 1), 1)
        call DAXPY(ilen, dummy, tmp4(srt), 1, hq(srt, 2), 1)
        call DAXPY(ilen, dummy, tmp3(srt), 1, hq(srt, 3), 1)
#else
        do ij = srt, end
            hq(ij, 1) = hq(ij, 1) - tmp2(ij)
            hq(ij, 2) = hq(ij, 2) - tmp4(ij)
            hq(ij, 3) = hq(ij, 3) - tmp3(ij)
        end do
#endif
!$omp end parallel
    end if

    ! #######################################################################
    ! Boundary conditions
    ! #######################################################################
    BcsFlowKmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowKmax%ref = 0.0_wp ! Scalar BCs initialized at start of routine

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Z(ibc, imax, jmax, kmax, g(3), hq(:, iq), &
                                        BcsFlowKmin%ref(:, :, iq), BcsFlowKmax%ref(:, :, iq), tmp1)
        end if

        p_hq(:, :, 1, iq) = BcsFlowKmin%ref(:, :, iq)
        p_hq(:, :, kmax, iq) = BcsFlowKmax%ref(:, :, iq)

    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Z(ibc, imax, jmax, kmax, g(3), hs(:, is), &
                                        BcsScalKmin%ref(:, :, is), BcsScalKmax%ref(:, :, is), tmp1)
        end if

        if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
            BcsScalKmax%type(is) /= DNS_SFC_STATIC) then
            call BOUNDARY_BCS_SURFACE_Z(is, s, hs, tmp1, tmp2)
        end if

        p_hs(:, :, 1, is) = BcsScalKmin%ref(:, :, is)
        p_hs(:, :, kmax, is) = BcsScalKmax%ref(:, :, is)

    end do

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_1
