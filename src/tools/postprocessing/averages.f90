#include "tlab_error.h"

program AVERAGES
    use TLab_Constants, only: wp, wi, small_wp, MAX_AVG_TEMPORAL
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_txc, isize_wrk3d, inb_flow, inb_scal, isize_field
    use TLab_Time, only: itime, rtime
    use TLab_Arrays
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Pointers, only: pointers_dt, u, v, w
#ifdef USE_MPI
    ! use mpi_f08, only: MPI_COMM_WORLD, MPI_REAL4
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: g, FDM_Initialize
    use NavierStokes
    use Thermodynamics, only: Thermo_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize, gravityProps, Gravity_Source, bbackground
    use SpecialForcing, only: SpecialForcing_Initialize
    ! use Rotation, only: Rotation_Initialize
    use Microphysics, only: Microphysics_Initialize
    use Radiation, only: Radiation_Initialize
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use OPR_Partial
    use OPR_Fourier
    use OPR_Elliptic
    use NSE_Burgers, only: NSE_Burgers_Initialize
    use NSE_Pressure
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN
    use Tensor
    use Reductions, only: Reduce_Block_InPlace

    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! iterations to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 30    ! options to be processed
    integer(wi) iopt_size
    integer(wi) opt_vec(iopt_size_max)
    ! real(wp) opt_vec2(iopt_size_max)

    integer, parameter :: params_size_max = 2

    ! -------------------------------------------------------------------
    ! Additional local arrays
    real(wp), allocatable, save :: mean(:), z_aux(:)

    integer(1), allocatable, save :: gate(:)
    type(pointers_dt) :: vars(16)
    ! real(wp), allocatable, save :: surface(:, :, :)       ! Gate envelopes

    character*32 fname !, varname(16)
    character*64 str

    integer opt_main, opt_block, opt_order
    integer nfield, ifield, is, ij, kmax_aux
    real(wp) eloc1, eloc2, eloc3, cos1, cos2, cos3, dummy
    logical iread_flow, iread_scal
    real(wp) params(2)
    ! integer(wi) io_sizes(5)
    ! type(io_subarray_dt) io_envelopes

    ! ! Gates for the definition of the intermittency function (partition of the fields)
    ! integer opt_cond, opt_cond_scal, opt_cond_relative
    ! integer, parameter :: igate_size_max = 8
    ! integer(wi) igate_size
    ! real(wp) gate_threshold(igate_size_max)
    ! integer(1) gate_level

    !########################################################################
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

    call Gravity_Initialize(ifile)
    call SpecialForcing_Initialize(ifile)
    ! call Rotation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)

    call TLab_Consistency_Check()

    call Averages_Initialize()

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call OPR_Partial_Initialize()
    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    call TLab_Initialize_Background(ifile)  ! Initialize thermodynamic quantities
    call NSE_Burgers_Initialize(ifile)

    ! do ig = 1, 3
    !     call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    ! end do

    allocate (gate(isize_field))

    ! in case g(2)%size is not divisible by opt_block, drop the upper most planes
    kmax_aux = z%size/opt_block
    allocate (z_aux(kmax_aux))          ! Reduced vertical grid
    z_aux(:) = 0.0_wp
    do ij = 1, kmax_aux*opt_block
        is = (ij - 1)/opt_block + 1
        z_aux(is) = z_aux(is) + z%nodes(ij)/real(opt_block, wp)
    end do

    if (opt_main == 1) then
        allocate (mean(kmax*MAX_AVG_TEMPORAL))
        ! else if (opt_main == 2) then
        !     allocate (mean(igate_size*(kmax_aux + 1)))
    else
        allocate (mean(opt_order*nfield*(kmax_aux + 1)))
    end if

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (str, *) itime; str = 'Processing iteration It'//trim(adjustl(str))
        call TLab_Write_ASCII(lfile, str)

        if (iread_scal) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
        end if

        if (iread_flow) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        call TLab_Diagnostic(imax, jmax, kmax, s)  ! Initialize diagnostic thermodynamic quantities

        ! ! -------------------------------------------------------------------
        ! ! Calculate intermittency
        ! ! -------------------------------------------------------------------
        ! if (opt_cond == 1) then ! External file
        !     write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname))
        !     call IO_Read_Field_INT1(fname, imax, jmax, kmax, itime, gate, params(1:2))
        !     igate_size = int(params(2))

        !     if (opt_main == 2) rtime = params(1)

        ! else if (opt_cond > 1) then
        !     opt_cond_scal = 1
        !     if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
        !         opt_cond_scal = inb_scal_array
        !     end if

        !     call TLab_Write_ASCII(lfile, 'Calculating partition...')
        !     call FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
        !                  imax, jmax, kmax, igate_size, gate_threshold, q, s, txc, gate)

        !     if (kmax_aux*opt_block /= g(2)%size) then
        !         call REDUCE_BLOCK_INPLACE_INT1(imax, jmax, kmax, i1, i1, i1, imax, kmax_aux*opt_block, kmax, gate, wrk1d)
        !     end if

        ! end if

        ! -------------------------------------------------------------------
        ! Type of averages
        ! -------------------------------------------------------------------
        select case (opt_main)

            ! ###################################################################
            ! Conventional statistics
            ! ###################################################################
        case (1)
            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                call NSE_Pressure_Incompressible(q, s, txc(:, 9), txc(:, 2), txc(:, 5), txc(:, 6))
            end if

            if (scal_on) then
                do is = 1, inb_scal_array          ! All, prognostic and diagnostic fields in array s
                    txc(1:isize_field, 6) = txc(1:isize_field, 9) ! Pass the pressure in tmp6
                    call AVG_SCAL_XZ(is, q, s, s(1, is), &
                                     txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), mean)
                end do

            end if

            if (flow_on) then
                txc(1:isize_field, 3) = txc(1:isize_field, 9) ! Pass the pressure in tmp3
                call AVG_FLOW_XZ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), &
                                 txc(1, 7), txc(1, 8), txc(1, 9), mean)
            end if

            ! ###################################################################
            ! Partition of field
            ! ###################################################################
        case (2)
            write (fname, *) itime; fname = 'int'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')

            ! do is = 1, igate_size
            !     write (varname(is), *) is; varname(is) = 'Partition'//trim(adjustl(varname(is)))
            ! end do
            ! call INTER_N_XZ(fname, itime, rtime, imax, jmax, kmax, igate_size, varname, gate, z%nodes, mean)

            ! if (opt_cond > 1) then ! write only if the gate information has not been read
            !     write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname))
            !     params(1) = rtime; params(2) = real(igate_size, wp)
            !     call IO_Write_Field_INT1(fname, imax, jmax, kmax, itime, gate, params(1:2))

            !     do is = 1, igate_size
            !         gate_level = int(is, KIND=1)
            !         call BOUNDARY_LOWER_INT1(imax, jmax, kmax, gate_level, z%nodes, gate, wrk3d, wrk2d, wrk2d(1, 2))
            !         do k = 1, kmax ! rearranging
            !             ij = (k - 1)*imax + 1
            !             surface(1:imax, is, k) = wrk2d(ij:ij + imax - 1, 1)
            !         end do
            !         call BOUNDARY_UPPER_INT1(imax, jmax, kmax, gate_level, z%nodes, gate, wrk3d, wrk2d, wrk2d(1, 2))
            !         do k = 1, kmax ! rearranging
            !             ij = (k - 1)*imax + 1
            !             surface(1:imax, is + igate_size, k) = wrk2d(ij:ij + imax - 1, 1)
            !         end do
            !     end do
            !     varname = ''
            !     write (fname, *) itime; fname = 'envelopesJ.'//trim(adjustl(fname))
            !     call IO_Write_Subarray(io_envelopes, fname, varname, surface, io_sizes)

            ! end if

            ! ###################################################################
            ! Momentum equation
            ! ###################################################################
        case (3)    ! To be checked
            write (fname, *) itime; fname = 'avgMom'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'U'
            ifield = ifield + 1; vars(ifield)%field => q(:, 2); vars(ifield)%tag = 'V'

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Uz'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Uzz'
            call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, g(3), q(1, 1), txc(1, 2), txc(1, 1))
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Vz'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Vzz'
            call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, g(3), q(1, 2), txc(1, 4), txc(1, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'WU)z'
            txc(1:isize_field, 6) = q(1:isize_field, 3)*q(1:isize_field, 1)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 6), txc(1, 5))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'WUz'
            txc(1:isize_field, 6) = q(1:isize_field, 3)*txc(1:isize_field, 1)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'UUx'
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), q(1, 1), txc(1, 7))
            txc(1:isize_field, 7) = q(1:isize_field, 1)*txc(1:isize_field, 7)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'VUy'
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), q(1, 1), txc(1, 8))
            txc(1:isize_field, 8) = q(1:isize_field, 2)*txc(1:isize_field, 8)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'WV)z'
            txc(1:isize_field, 10) = q(1:isize_field, 2)*q(1:isize_field, 3)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 10), txc(1, 9))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 10); vars(ifield)%tag = 'WVz'
            txc(1:isize_field, 10) = q(1:isize_field, 3)*txc(1:isize_field, 3)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 11); vars(ifield)%tag = 'UVx'
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), q(1, 2), txc(1, 11))
            txc(1:isize_field, 11) = q(1:isize_field, 1)*txc(1:isize_field, 11)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 12); vars(ifield)%tag = 'VVy'
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), q(1, 2), txc(1, 12))
            txc(1:isize_field, 12) = q(1:isize_field, 3)*txc(1:isize_field, 12)

            ! ###################################################################
            ! Main variables
            ! ###################################################################
        case (4)
            write (fname, *) itime; fname = 'avgMain'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ifield = ifield + 1; vars(ifield)%field => u(:); vars(ifield)%tag = 'U'
            ifield = ifield + 1; vars(ifield)%field => v(:); vars(ifield)%tag = 'V'
            ifield = ifield + 1; vars(ifield)%field => w(:); vars(ifield)%tag = 'W'

            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
                ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'P'
            else
                ifield = ifield + 1; vars(ifield)%field => q(:, 5); vars(ifield)%tag = 'R'
                ifield = ifield + 1; vars(ifield)%field => q(:, 6); vars(ifield)%tag = 'P'
                ifield = ifield + 1; vars(ifield)%field => q(:, 7); vars(ifield)%tag = 'T'
            end if

            if (scal_on) then
                do is = 1, inb_scal_array          ! All, prognostic and diagnostic fields in array s
                    ifield = ifield + 1; vars(ifield)%field => s(:, is); write (vars(ifield)%tag, *) is; vars(ifield)%tag = 'Scalar'//trim(adjustl(vars(ifield)%tag))
                end do
            end if

            ! ###################################################################
            ! Enstrophy equation
            ! ###################################################################
        case (5)
            write (fname, *) itime; fname = 'avgW2'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ! result vector in txc4, txc5, txc6
            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                ! txc(:, 4) = 0.0_wp; txc(:, 5) = 0.0_wp; txc(:, 6) = 0.0_wp

                select case (nse_eqns)
                case (DNS_EQNS_ANELASTIC)
                    ! call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, wrk3d)

                case (DNS_EQNS_BOUSSINESQ)
                    wrk1d(1:kmax, 1) = bbackground(1:kmax)
                    bbackground(1:kmax) = 0.0_wp
                    call Gravity_Source(gravityProps, imax, jmax, kmax, s, wrk3d)
                    bbackground(1:kmax) = wrk1d(1:kmax, 1)
                end select
                s(1:isize_field, 1) = wrk3d(1:isize_field)*gravityProps%vector(2)

                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s, txc(1, 4))
                txc(:, 4) = -txc(:, 4)
                txc(:, 5) = 0.0_wp
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s, txc(1, 6))

            else
                call FI_VORTICITY_BAROCLINIC(imax, jmax, kmax, q(1, 5), q(1, 6), txc(1, 4), txc(1, 3), txc(1, 7))
            end if

            call FI_CURL(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 7))
            txc(1:isize_field, 8) = txc(1:isize_field, 1)*txc(1:isize_field, 4) &
                                    + txc(1:isize_field, 2)*txc(1:isize_field, 5) + txc(1:isize_field, 3)*txc(1:isize_field, 6)

            call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, u, v, w, txc(1, 1), &
                                         txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

            call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, u, v, w, txc(1, 2), &
                                        txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = visc*txc(1:isize_field, 2)

            call FI_VORTICITY(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4), txc(1, 5))  ! Enstrophy
            call FI_INVARIANT_P(imax, jmax, kmax, u, v, w, txc(1, 4), txc(1, 5))  ! Dilatation

            txc(1:isize_field, 6) = txc(1:isize_field, 4)*txc(1:isize_field, 3) ! -w^2 div(u)
            txc(1:isize_field, 5) = txc(1:isize_field, 1)/txc(1:isize_field, 3) ! production rate
            txc(1:isize_field, 4) = log(txc(1:isize_field, 3))                  ! ln(w^2)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'EnstrophyW_iW_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'LnEnstrophyW_iW_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionW_iW_jS_ij'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuW_iLapW_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'DilatationMsW_iW_iDivU'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Baroclinic'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'RateAN_iN_jS_ij'

            ! ###################################################################
            ! Strain equation
            ! ###################################################################
        case (6)
            write (fname, *) itime; fname = 'avgS2'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, txc(1, 1), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            else
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, q(1, 6), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            end if
            txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 2)

            call FI_STRAIN_PRODUCTION(imax, jmax, kmax, u, v, w, &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = 2.0_wp*txc(1:isize_field, 2)

            call FI_STRAIN_DIFFUSION(imax, jmax, kmax, u, v, w, &
                                     txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8))
            txc(1:isize_field, 3) = 2.0_wp*visc*txc(1:isize_field, 3)

            call FI_STRAIN(imax, jmax, kmax, u, v, w, txc(1, 4), txc(1, 5), txc(1, 6))
            txc(1:isize_field, 4) = 2.0_wp*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = log(txc(1:isize_field, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Strain2S_ijS_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'LnStrain2S_ijS_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'ProductionMs2S_ijS_jkS_ki'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'DiffusionNuS_ijLapS_ij'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Pressure2S_ijP_ij'

            ! ###################################################################
            ! Scalar gradient equation
            ! ###################################################################
        case (7)
            write (fname, *) itime; fname = 'avgG2'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s, u, v, w, &
                                        txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

            call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s, &   ! array u used as auxiliar
                                       txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), u)
            txc(1:isize_field, 2) = txc(1:isize_field, 2)*visc/schmidt(inb_scal)

            call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 3), txc(1, 4))
            txc(1:isize_field, 5) = txc(1:isize_field, 1)/txc(1:isize_field, 3)
            txc(1:isize_field, 4) = log(txc(1:isize_field, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'GradientG_iG_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'LnGradientG_iG_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionMsG_iG_jS_ij'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuG_iLapG_i'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'StrainAMsN_iN_jS_ij'

            ! ###################################################################
            ! Velocity gradient invariants
            ! ###################################################################
        case (8)
            write (fname, *) itime; fname = 'avgInv'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_INVARIANT_R(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call FI_INVARIANT_Q(imax, jmax, kmax, u, v, w, txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
            call FI_INVARIANT_P(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'InvariantP'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'InvariantQ'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'InvariantR'

            ! ###################################################################
            ! Scalar gradient components
            ! ###################################################################
        case (9)
            write (fname, *) itime; fname = 'avgGi'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s, txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s, txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s, txc(1, 3))
            do ij = 1, isize_field                       ! Angles; s array is overwritten to save space
                dummy = txc(ij, 3)/sqrt(txc(ij, 1)*txc(ij, 1) + txc(ij, 2)*txc(ij, 2) + txc(ij, 3)*txc(ij, 3))
                txc(ij, 4) = asin(dummy)                  ! with Oz
                s(ij, 1) = atan2(txc(ij, 2), txc(ij, 1))    ! with Ox in plane xOy
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'GradientX'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'GradientY'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'GradientZ'
            ifield = ifield + 1; vars(ifield)%field => s(:, 1); vars(ifield)%tag = 'Theta'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Phi'

            ! ###################################################################
            ! eigenvalues of rate-of-strain tensor
            ! ###################################################################
        case (10)
            write (fname, *) itime; fname = 'avgEig'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_STRAIN_TENSOR(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'Lambda1'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Lambda2'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'Lambda3'

            ! ###################################################################
            ! eigenframe of rate-of-strain tensor
            ! ###################################################################
        case (11)
            write (fname, *) itime; fname = 'avgCos'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_STRAIN_TENSOR(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7))  ! txc7-txc9
            call TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1, 1), txc(1, 7))   ! txc1-txc6

            call FI_CURL(imax, jmax, kmax, u, v, w, txc(1, 7), txc(1, 8), txc(1, 9), txc(1, 10))
            do ij = 1, isize_field                                             ! local direction cosines of vorticity vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                u(ij) = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                v(ij) = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                w(ij) = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
            end do

            ifield = ifield + 1; vars(ifield)%field => u; vars(ifield)%tag = 'cos(w,lambda1)'
            ifield = ifield + 1; vars(ifield)%field => v; vars(ifield)%tag = 'cos(w,lambda2)'
            ifield = ifield + 1; vars(ifield)%field => w; vars(ifield)%tag = 'cos(w,lambda3)'

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s, txc(1, 7))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s, txc(1, 8))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s, txc(1, 9))
            do ij = 1, isize_field                                             ! local direction cosines of scalar gradient vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                cos1 = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                cos2 = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                cos3 = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
                txc(ij, 7) = cos1; txc(ij, 8) = cos2; txc(ij, 9) = cos3
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'cos(G,lambda1)'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'cos(G,lambda2)'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'cos(G,lambda3)'

            ! ###################################################################
            ! longitudinal velocity derivatives
            ! ###################################################################
        case (12)
            write (fname, *) itime; fname = 'avgDer'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), u, txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), v, txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), w, txc(1, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'dudx'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'dvdy'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'dwdz'

            ! ###################################################################
            ! Vertical fluxes
            ! ###################################################################
        case (13)
            ifield = 0
            write (fname, *) itime; fname = 'avgFluxZ'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), u, txc(:, 1))
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), v, txc(:, 2))
            txc(:, 1) = (txc(:, 1) + txc(:, 2))*visc
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'tauzx'

            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), w, txc(:, 2))
            txc(:, 2) = txc(:, 2)*2.0_wp*visc
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'tauzz'

            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), w, txc(:, 3))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), v, txc(:, 4))
            txc(:, 3) = (txc(:, 3) + txc(:, 4))*visc
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'tauzy'

            do is = 1, inb_scal_array
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s(:, is), txc(:, 3 + is))
                txc(:, 3 + is) = txc(:, 3 + is)*visc/schmidt(inb_scal)
                ifield = ifield + 1; vars(ifield)%field => txc(:, 3 + is); write (vars(ifield)%tag, *) is; vars(ifield)%tag = 'tauz'//trim(adjustl(vars(ifield)%tag))
            end do

            u = u*w
            ifield = ifield + 1; vars(ifield)%field => u; vars(ifield)%tag = 'wu'
            v = v*w
            ifield = ifield + 1; vars(ifield)%field => v; vars(ifield)%tag = 'wv'
            ! I need w below
            ifield = ifield + 1; vars(ifield)%field => w; vars(ifield)%tag = 'ww'
            do is = 1, inb_scal_array
                s(:, is) = s(:, is)*w
                ifield = ifield + 1; vars(ifield)%field => s(:, is); write (vars(ifield)%tag, *) is; vars(ifield)%tag = 'w'//trim(adjustl(vars(ifield)%tag))
            end do
            w = w*w ! I need w above for the scalar fluxes

            ! ###################################################################
            ! Hydrostatic and dynamic pressure
            ! ###################################################################
        case (14)
            write (fname, *) itime; fname = 'avgP'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'P'

            q = 0.0_wp
            call NSE_Pressure_Incompressible(q, s, txc(:, 2), txc(:, 3), txc(:, 6), txc(:, 7))
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Psta'

            txc(:, 3) = txc(:, 1) - txc(:, 2)
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Pdyn'

            ! ###################################################################
            ! Dissipation
            ! ###################################################################
        case (15)
            write (fname, *) itime; fname = 'avgEps'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            ! call FI_DISSIPATION(i1, imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
            ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*visc

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Eps'

            ! ###################################################################
            ! Covariances among scalars
            ! ###################################################################
        case (16)
            write (fname, *) itime; fname = 'avgSiCov'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, s(1, 1))
            call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, s(1, 2))

            txc(1:isize_field, 1) = s(1:isize_field, 1)*s(1:isize_field, 2)
            txc(1:isize_field, 2) = txc(1:isize_field, 1)*s(1:isize_field, 1)
            txc(1:isize_field, 3) = txc(1:isize_field, 1)*s(1:isize_field, 2)

            ifield = 0
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 's1s2'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 's1s2s1'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 's1s2s2'

            ! ###################################################################
            ! Potential vorticity
            ! ###################################################################
        case (17)
            write (fname, *) itime; fname = 'avgPV'//trim(adjustl(fname))
            call TLab_Write_ASCII(lfile, 'Computing '//trim(adjustl(fname))//'...')
            ifield = 0

            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
            txc(1:isize_field, 6) = txc(1:isize_field, 1)*txc(1:isize_field, 1) &
                                    + txc(1:isize_field, 2)*txc(1:isize_field, 2) &
                                    + txc(1:isize_field, 3)*txc(1:isize_field, 3) ! Enstrophy
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 3)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b

            txc(1:isize_field, 5) = sqrt(txc(1:isize_field, 5) + small_wp)
            txc(1:isize_field, 6) = sqrt(txc(1:isize_field, 6) + small_wp)
            txc(1:isize_field, 2) = txc(1:isize_field, 1)/(txc(1:isize_field, 5)*txc(1:isize_field, 6)) ! Cosine of angle between 2 vectors

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'PV'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Cos'

        end select

        if (opt_main > 2) then
            if (nfield < ifield) then ! Check
                call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
                call TLab_Stop(DNS_ERROR_WRKSIZE)
            end if

            if (kmax_aux*opt_block /= g(2)%size) then
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, kmax_aux*opt_block, kmax, vars(is)%field)
                end do
            end if

            ! call AVG_N_XZ(fname, itime, rtime, imax*opt_block, kmax_aux, kmax, &
            !               ifield, opt_order, vars, gate_level, gate, z_aux, mean)

        end if

    end do
    call TLab_Stop(0)

contains
    subroutine Averages_Initialize()

        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        ! #######################################################################
        ! Read from tlab.ini
        bakfile = trim(adjustl(ifile))//'.bak'
        block = 'PostProcessing'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        ! -------------------------------------------------------------------
        call ScanFile_Char(bakfile, ifile, block, 'Files', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Iteration numbers ?'
            read (*, '(A512)') sRes
#endif
        end if
        itime_size = itime_size_max
        call LIST_INTEGER(sRes, itime_size, itime_vec)

        if (itime_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing Files input.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! -------------------------------------------------------------------
        opt_main = -1 ! default values
        opt_block = 1
        ! gate_level = 0
        opt_order = 1

        call ScanFile_Char(bakfile, ifile, block, 'ParamAverages', '-1', sRes)
        iopt_size = iopt_size_max
        call LIST_INTEGER(sRes, iopt_size, opt_vec)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Option ?'
            write (*, *) ' 1. Conventional averages'
            write (*, *) ' 2. Intermittency or gate function'
            write (*, *) ' 3. Momentum equation'
            write (*, *) ' 4. Main variables'
            write (*, *) ' 5. Enstrophy W_iW_i/2 equation'
            write (*, *) ' 6. Strain 2S_ijS_ij/2 equation'
            write (*, *) ' 7. Scalar gradient G_iG_i/2 equation'
            write (*, *) ' 8. Velocity gradient invariants'
            write (*, *) ' 9. Scalar gradient components'
            write (*, *) '10. Eigenvalues of rate-of-strain tensor'
            write (*, *) '11. Eigenframe of rate-of-strain tensor'
            write (*, *) '12. Longitudinal velocity derivatives'
            write (*, *) '13. Vertical fluxes'
            write (*, *) '14. Pressure partition'
            write (*, *) '15. Dissipation'
            write (*, *) '16. Third-order scalar covariances'
            write (*, *) '17. Potential vorticity'
            read (*, *) opt_main

            write (*, *) 'Planes block size ?'
            read (*, *) opt_block

            if (opt_main > 2) then
                ! write (*, *) 'Gate level to be used ?'
                ! read (*, *) gate_level
                write (*, *) 'Number of moments ?'
                read (*, *) opt_order
            end if

#endif
        else
            opt_main = int(opt_vec(1))
            if (iopt_size >= 2) opt_block = int(opt_vec(2))
            if (iopt_size >= 3) opt_order = int(opt_vec(4))
            ! if (iopt_size >= 3) gate_level = int(opt_vec(3), KIND=1)

        end if

        if (opt_main < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing input [ParamAverages] in tlab.ini.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        if (opt_block < 1) then
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Invalid value of opt_block.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

! -------------------------------------------------------------------
! Defining gate levels for conditioning
! -------------------------------------------------------------------
!         opt_cond = 0 ! default values
!         opt_cond_relative = 0
        ! igate_size = 0

!         if (opt_main == 2 .or. gate_level /= 0) then
! #include "dns_read_partition.h"
!             if (opt_cond > 1) inb_txc = max(inb_txc, 5)
!         end if

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        iread_flow = .false.
        iread_scal = .false.
        inb_txc = 0
        nfield = 2

        select case (opt_main)
        case (1)
            inb_txc = max(inb_txc, 9)
            iread_flow = flow_on; iread_scal = scal_on
        case (2)
        case (3)
            nfield = 14
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 12)
        case (4)
            nfield = 6 + inb_scal
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 3)
            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) inb_txc = max(inb_txc, 6)
        case (5) ! enstrophy
            nfield = 7
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 8)
        case (6)
            nfield = 5
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 8)
        case (7) ! scalar gradient
            nfield = 5
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 6)
        case (8)
            nfield = 3
            iread_flow = .true.; inb_txc = max(inb_txc, 6)
        case (9)
            nfield = 5
            iread_scal = .true.; inb_txc = max(inb_txc, 4)
        case (10) ! eigenvalues
            nfield = 3
            iread_flow = .true.; inb_txc = max(inb_txc, 9)
        case (11) ! eigenframe
            nfield = 6
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 10)
        case (12) ! longitudinal velocity derivatives
            nfield = 3
            iread_flow = .true.; iread_scal = .false.; inb_txc = max(inb_txc, 3)
        case (13) ! Vertical flux
            nfield = 2*(3 + inb_scal_array)
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(max(inb_txc, 3 + inb_scal_array), 4)
        case (14) ! pressure partition
            nfield = 3
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7)
        case (15) ! dissipation partition
            nfield = 1
            iread_flow = .true.; iread_scal = .false.; inb_txc = max(inb_txc, 6)
        case (16) ! third-order scalar covariances
            nfield = 3
            iread_flow = .false.; iread_scal = .true.; inb_txc = max(inb_txc, 3)
        case (17) ! potential vorticity
            nfield = 2
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 6)
        end select

        isize_wrk3d = max(isize_wrk3d, opt_order*nfield*jmax)

        ! #######################################################################
        ! Subarray information to read and write envelopes
!         if (opt_main == 2 .and. opt_cond > 1 .and. igate_size /= 0) then
!             ! Info for IO routines: total size, lower bound, upper bound, stride, # variables
!             idummy = imax*2*igate_size*kmax; io_sizes = [idummy, 1, idummy, 1, 1]
!             io_envelopes%offset = 0
!             io_envelopes%precision = IO_TYPE_SINGLE
! #ifdef USE_MPI
!             io_envelopes%active = .true.
!             io_envelopes%communicator = MPI_COMM_WORLD
!             io_envelopes%subarray = IO_Create_Subarray_XOZ(imax, igate_size*2, kmax, MPI_REAL4)
! #endif

!             allocate (surface(imax, 2*igate_size, kmax))

!         end if

        return
    end subroutine Averages_Initialize

end program AVERAGES
