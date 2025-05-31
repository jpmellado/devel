! Calculating pressure in incompressible modes
! Similar to corresponding routines in time marching module, but reorganize
! to save memory at the expense of time, since these routines are only diagnostic

module NSE_Pressure
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_txc_field
    use TLab_Memory, only: inb_flow_array, inb_flow, inb_scal_array
    use FDM, only: g
    use OPR_Partial
    use OPR_Elliptic
    use NavierStokes
    use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_InPlace
    use NSE_Burgers
    use TLab_Sources
    implicit none
    private

    public :: NSE_Pressure_Incompressible

    real(wp), allocatable :: bcs_hb(:, :), bcs_ht(:, :)
    real(wp), pointer :: u(:), v(:), w(:), p_hq3(:, :, :)

contains
    ! #######################################################################
    ! #######################################################################
    subroutine NSE_Pressure_Incompressible(q, s, p, hq, tmp1, tmp2)
        real(wp), intent(in), target :: q(isize_field, inb_flow_array)
        real(wp), intent(in) :: s(isize_field, inb_scal_array)
        real(wp), intent(out) :: p(isize_field)
        real(wp), intent(inout), target :: hq(isize_field, inb_flow)
        real(wp), intent(inout) :: tmp1(isize_txc_field), tmp2(isize_txc_field)

        ! #######################################################################
        ! Define pointers
        u => q(:, 1)
        v => q(:, 2)
        w => q(:, 3)

        p_hq3(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, 3)
        
        ! #######################################################################
        ! Diffusion and advection terms
        ! Using p as auxiliary array
        ! #######################################################################
        call NSE_Burgers_X(0, imax, jmax, kmax, u, hq(:, 1), p)          ! store u transposed in p
        call NSE_Burgers_X(0, imax, jmax, kmax, v, hq(:, 2), tmp2, p)    ! p contains u transposed
        call NSE_Burgers_X(0, imax, jmax, kmax, w, hq(:, 3), tmp2, p)    ! p contains u transposed

        call NSE_Burgers_Y(0, imax, jmax, kmax, u, tmp1, p)              ! store v transposed in p
        hq(:, 1) = hq(:, 1) + tmp1(:)
        call NSE_Burgers_Y(0, imax, jmax, kmax, v, tmp1, tmp2, p)        ! p contains v transposed
        hq(:, 2) = hq(:, 2) + tmp1(:)
        call NSE_Burgers_Y(0, imax, jmax, kmax, w, tmp1, tmp2, p)        ! p contains v transposed
        hq(:, 3) = hq(:, 3) + tmp1(:)

        call NSE_Burgers_Z(0, imax, jmax, kmax, u, tmp1, w)
        hq(:, 1) = hq(:, 1) + tmp1(:)
        call NSE_Burgers_Z(0, imax, jmax, kmax, v, tmp1, w)
        hq(:, 2) = hq(:, 2) + tmp1(:)
        call NSE_Burgers_Z(0, imax, jmax, kmax, w, tmp1, w)
        hq(:, 3) = hq(:, 3) + tmp1(:)

        ! #######################################################################
        ! Add forcing terms
        ! #######################################################################
        call TLab_Sources_Flow(q, s, hq, tmp1)

        ! We are missing the buffer nudging here.

        ! #######################################################################
        ! Pressure term
        ! #######################################################################
        ! Forcing term
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, hq(:, 1))
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, hq(:, 2))
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, hq(:, 3))
        end if
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), hq(:, 1), p)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), hq(:, 2), tmp1)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), hq(:, 3), tmp2)
        p(:) = p(:) + tmp1(:) + tmp2(:)

        ! Neumman BCs in d/dz(p) s.t. v=0 (no-penetration)
        if (allocated(bcs_hb)) deallocate (bcs_hb)
        allocate (bcs_hb(imax, jmax))
        if (allocated(bcs_ht)) deallocate (bcs_ht)
        allocate (bcs_ht(imax, jmax))

        bcs_hb(:, :) = p_hq3(:, :, 1)
        bcs_ht(:, :) = p_hq3(:, :, kmax)

        call OPR_Poisson(imax, jmax, kmax, BCS_NN, p, tmp1, tmp2, bcs_hb, bcs_ht)

        return
    end subroutine NSE_Pressure_Incompressible

end module NSE_Pressure
