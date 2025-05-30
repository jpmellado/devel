!########################################################################
!# Sources (processes) in the evolution equations.
!########################################################################
module TLab_Sources
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_scal, inb_scal_array
    use TLab_OpenMP
    use FDM, only: g
    ! use FDM, only: fdm_Int0
    use NavierStokes, only: nse_eqns, DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC
    use Thermo_Anelastic, only: ribackground, Thermo_Anelastic_Buoyancy, Thermo_Anelastic_WEIGHT_ADD
    use Gravity, only: gravityProps, Gravity_Source
    ! use Rotation, only: coriolis, Rotation_Coriolis
    ! use Radiation
    use Microphysics
    use SpecialForcing
    ! use LargeScaleForcing
    implicit none
    private

    public :: TLab_Sources_Flow
    public :: TLab_Sources_Scal

contains
! #######################################################################
! #######################################################################
    subroutine TLab_Sources_Flow(q, s, hq, tmp1)
        use TLab_Time, only: rtime
        real(wp), intent(in) :: q(isize_field, *), s(isize_field, *)
        real(wp), intent(out) :: hq(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field)

        ! -----------------------------------------------------------------------
        integer iq

        ! #######################################################################
        ! call Rotation_Coriolis(coriolis, imax, jmax, kmax, q, hq)

        do iq = 1, 3
            ! -----------------------------------------------------------------------
            if (gravityProps%active(iq)) then
                select case (nse_eqns)
                case(DNS_EQNS_BOUSSINESQ)
                    call Gravity_Source(gravityProps, imax, jmax, kmax, s, tmp1)

                case(DNS_EQNS_ANELASTIC)
                    call Thermo_Anelastic_Buoyancy(imax, jmax, kmax, s, tmp1)

                end select
                
                hq(:, iq) = hq(:, iq) + gravityProps%vector(iq)*tmp1(:)

            end if

!             ! -----------------------------------------------------------------------
!             ! Subsidence
!             ! -----------------------------------------------------------------------
!             if (subsidenceProps%active(iq)) then
!                 call LargeScaleForcing_Subsidence(subsidenceProps, imax, jmax, kmax, q(:, iq), tmp1)

! !$omp parallel default( shared ) &
! !$omp private( ij, srt,end,siz )
!                 call TLab_OMP_PARTITION(isize_field, srt, end, siz)

!                 do ij = srt, end
!                     hq(ij, iq) = hq(ij, iq) + tmp1(ij)
!                 end do
! !$omp end parallel

!             end if

            ! -----------------------------------------------------------------------
            if (forcingProps%active(iq)) then
                call SpecialForcing_Source(forcingProps, imax, jmax, kmax, iq, rtime, q(:, iq), hq(:, iq), tmp1)

                hq(:, iq) = hq(:, iq) + tmp1(:)

            end if

        end do

        return
    end subroutine TLab_Sources_Flow

    ! #######################################################################
    ! #######################################################################
    subroutine TLab_Sources_Scal(s, hs, tmp1, tmp2, tmp3, tmp4)
        real(wp), intent(in) :: s(isize_field, *)
        real(wp), intent(out) :: hs(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field), tmp2(isize_field), tmp3(isize_field), tmp4(isize_field)

        ! -----------------------------------------------------------------------
        integer is

        ! #######################################################################
        do is = 1, inb_scal ! Start loop over the N scalars

!             ! -----------------------------------------------------------------------
!             ! Radiation
!             ! -----------------------------------------------------------------------
!             if (infraredProps%active(is)) then
!                 call Radiation_Infrared_Z(infraredProps, imax, jmax, kmax, fdm_Int0, s, tmp1, tmp2, tmp3, tmp4)

!                 if (nse_eqns == DNS_EQNS_ANELASTIC) then
!                     call Thermo_Anelastic_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(:, is))
!                 else
! !$omp parallel default( shared ) &
! !$omp private( ij, srt,end,siz )
!                     call TLab_OMP_PARTITION(isize_field, srt, end, siz)

!                     do ij = srt, end
!                         hs(ij, is) = hs(ij, is) + tmp1(ij)
!                     end do
! !$omp end parallel
!                 end if

!             end if

            ! -----------------------------------------------------------------------
            ! Microphysics
            ! -----------------------------------------------------------------------
            if (sedimentationProps%active(is)) then
                call Microphysics_Sedimentation(sedimentationProps, imax, jmax, kmax, is, g(3), s, tmp1, tmp2)

                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(:, is))
                else
                    hs(:, is) = hs(:, is) + tmp1(:)
                end if

            end if

!             ! -----------------------------------------------------------------------
!             ! Chemistry
!             ! -----------------------------------------------------------------------
!             if (chemistryProps%active(is)) then
!                 call Chemistry_Source(chemistryProps, imax, jmax, kmax, is, s, tmp1)

! !$omp parallel default( shared ) &
! !$omp private( ij, srt,end,siz )
!                 call TLab_OMP_PARTITION(isize_field, srt, end, siz)

!                 do ij = srt, end
!                     hs(ij, is) = hs(ij, is) + tmp1(ij)
!                 end do
! !$omp end parallel

!             end if

!             ! -----------------------------------------------------------------------
!             ! Subsidence
!             ! -----------------------------------------------------------------------
!             if (subsidenceProps%active(is)) then
!                 call LargeScaleForcing_Subsidence(subsidenceProps, imax, jmax, kmax, s(:, is), tmp1)

! !$omp parallel default( shared ) &
! !$omp private( ij, srt,end,siz )
!                 call TLab_OMP_PARTITION(isize_field, srt, end, siz)

!                 do ij = srt, end
!                     hs(ij, is) = hs(ij, is) + tmp1(ij)
!                 end do
! !$omp end parallel

!             end if

        end do

        return
    end subroutine TLab_Sources_Scal

end module TLab_Sources
