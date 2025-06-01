#include "tlab_error.h"

program VISUALS
    use TLab_Constants, only: wp, wi, small_wp, MAX_PARS, fmt_r
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_txc, inb_flow, inb_scal, isize_field
    use TLab_Arrays
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use mpi_f08, only: MPI_COMM_WORLD, MPI_REAL4
    use TLabMPI_VARS, only: ims_pro_i, ims_pro_j, ims_comm_x, ims_comm_y
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: g, FDM_Initialize
    use FDM, only: fdm_Int0
    use NavierStokes!, only: NavierStokes_Initialize_Parameters
    use Thermodynamics, only: Thermo_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize, gravityProps, Gravity_Source, bbackground
    ! use Rotation, only: Rotation_Initialize
    use Thermo_Anelastic
    use Radiation !, only: Radiation_Initialize, infraredProps
    use Microphysics !, only: Microphysics_Initialize, sedimentationProps
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

    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! iterations to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 30    ! options to be processed
    integer(wi) iopt_size, iv
    integer(wi) opt_vec(iopt_size_max)
    character(len=64) opt_name(iopt_size_max)

    integer :: opt_format                           ! File format
    integer, parameter :: FORMAT_SINGLE = 1         ! Single precision, no headers
    integer, parameter :: FORMAT_GENERAL = 2        ! General IO format

    integer(wi) subdomain(6)                        ! Subdomain to be saved

    character(len=32) time_str                      ! Time stamp
    integer, parameter :: MaskSize = 6

    character*32 flow_file, scal_file, plot_file
    character*64 str

    integer, parameter :: IO_SUBARRAY_VISUALS_XOY = 1
    integer, parameter :: IO_SUBARRAY_VISUALS_YOZ = 2
    integer, parameter :: IO_SUBARRAY_VISUALS_XOZ = 3
    type(io_subarray_dt) :: io_subarrays(3)

    integer(wi) ij, is
    integer(wi), parameter :: iscal_offset = 9 ! to be removed
    logical iread_flow, iread_scal
    real(wp) diff
    real(wp) params(MAX_PARS)

    ! ! Gates for the definition of the intermittency function (partition of the fields)
    ! integer(wi) opt_cond, opt_cond_scal, opt_cond_relative
    ! integer(wi), parameter :: igate_size_max = 8
    ! integer(wi) igate_size, ig
    ! real(wp) gate_threshold(igate_size_max)
    ! integer(1), dimension(:), allocatable, save :: gate

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
    ! call Rotation_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)

    call TLab_Consistency_Check()

    call Visuals_Initialize()

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    ! allocate (gate(isize_field))

    call OPR_Partial_Initialize()
    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    call TLab_Initialize_Background(ifile)
    call NSE_Burgers_Initialize(ifile)

    ! allocate (gate(isize_field))

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (str, *) itime; str = 'Processing iteration It'//trim(adjustl(str))//'.'
        call TLab_Write_ASCII(lfile, str)

        if (scal_on .and. iread_scal) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_Read_Fields(scal_file, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
        elseif (.not. scal_on) then
            s = 0.0_wp
        end if

        if (iread_flow) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_Read_Fields(flow_file, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        call TLab_Diagnostic(imax, jmax, kmax, s)

        write (str, fmt_r) rtime; str = 'Physical time '//trim(adjustl(str))
        call TLab_Write_ASCII(lfile, str)

        ! ! -------------------------------------------------------------------
        ! ! Calculate intermittency
        ! ! -------------------------------------------------------------------
        ! if (opt_cond == 1) then ! Read external file
        !     write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname))
        !     call IO_Read_Field_INT1(fname, imax, jmax, kmax, itime, gate, params(1:2))
        !     igate_size = int(params(2))

        ! else if (opt_cond > 1) then
        !     opt_cond_scal = 1 ! Scalar field to use for the conditioning
        !     if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
        !         opt_cond_scal = inb_scal_array
        !     end if

        !     call TLab_Write_ASCII(lfile, 'Calculating partition...')
        !     call FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
        !                  imax, jmax, kmax, igate_size, gate_threshold, q, s, txc, gate)
        ! end if

        ! -------------------------------------------------------------------
        ! define time string
        ! -------------------------------------------------------------------
        do ij = MaskSize, 1, -1
            time_str(ij:ij) = '0'
        end do
        write (plot_file, '(I10)') itime
        time_str(MaskSize - len_trim(adjustl(plot_file)) + 1:Masksize) = trim(adjustl(plot_file))

        ! -------------------------------------------------------------------
        ! Loop over options
        ! -------------------------------------------------------------------
        do iv = 1, iopt_size
            plot_file = trim(adjustl(opt_name(opt_vec(iv))))//time_str(1:MaskSize)

            select case (trim(adjustl(opt_name(opt_vec(iv)))))

                ! ###################################################################
                ! Velocities
                ! ###################################################################
            case ('VelocityX')
                txc(1:isize_field, 1) = q(1:isize_field, 1)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('VelocityY')
                txc(1:isize_field, 1) = q(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('VelocityZ')
                txc(1:isize_field, 1) = q(1:isize_field, 3)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('VelocityVector')
                txc(1:isize_field, 1:3) = q(1:isize_field, 1:3)
                call Write_Visuals(plot_file, txc(:, 1:3))

            case ('VelocityMagnitude')
                txc(1:isize_field, 1) = sqrt(Tensor_Dot(q(1:isize_field, 1:3), q(1:isize_field, 1:3)))
                ! sqrt(q(1:isize_field, 1)**2 + q(1:isize_field, 2)**2 + q(1:isize_field, 3)**2)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! ###################################################################
                ! Thermodynamic state
                ! ###################################################################
            case ('Density')
                select case (nse_eqns)
                case (DNS_EQNS_COMPRESSIBLE)
                    txc(1:isize_field, 1) = q(1:isize_field, 5)

                case (DNS_EQNS_ANELASTIC)
                    call Thermo_Anelastic_Rho(imax, jmax, kmax, s, txc(:, 1), wrk3d)

                case (DNS_EQNS_BOUSSINESQ)      ! Using buoyancy to calculate density
                    wrk1d(1:kmax, 1) = bbackground(1:kmax)
                    bbackground(1:kmax) = 0.0_wp
                    call Gravity_Source(gravityProps, imax, jmax, kmax, s, txc(1, 1))
                    txc(1:isize_field, 1) = txc(1:isize_field, 1)/froude + 1.0_wp
                    bbackground(1:kmax) = wrk1d(1:kmax, 1)

                end select
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('Temperature')
                select case (nse_eqns)
                case (DNS_EQNS_COMPRESSIBLE)
                    txc(1:isize_field, 1) = q(1:isize_field, 7)

                case (DNS_EQNS_ANELASTIC)
                    call Thermo_Anelastic_T(imax, jmax, kmax, s, txc(1, 1))

                end select
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('Pressure')
                select case (nse_eqns)
                case (DNS_EQNS_COMPRESSIBLE)
                    txc(1:isize_field, 2) = q(1:isize_field, 6)

                case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)      ! Although Pressure is not thermodynamic variable here...
                    call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
                    txc(1:isize_field, 2) = txc(1:isize_field, 1)   ! Save txc1 for later caculations

                end select
                call Write_Visuals(plot_file, txc(:, 2:2))

                ! Additional pressure stuff...
                plot_file = 'PressureGradientPower'//time_str(1:MaskSize)
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), txc(1, 1), txc(1, 2))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 1), txc(1, 3))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 1), txc(1, 4))
                txc(1:isize_field, 2) = -Tensor_Dot(q(:, 1:3), txc(:, 2:4))
                call Write_Visuals(plot_file, txc(:, 2:2))

                txc(1:isize_field, 2) = txc(1:isize_field, 1)
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(1, 2))

                plot_file = 'PressureStrainX'//time_str(1:MaskSize)
                txc(1:isize_field, 3) = q(1:isize_field, 1)
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), txc(1, 3), txc(1, 4))
                txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                call Write_Visuals(plot_file, txc(:, 3:3))

                plot_file = 'PressureStrainY'//time_str(1:MaskSize)
                txc(1:isize_field, 3) = q(1:isize_field, 2)
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 3), txc(1, 4))
                txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                call Write_Visuals(plot_file, txc(:, 3:3))

                plot_file = 'PressureStrainZ'//time_str(1:MaskSize)
                txc(1:isize_field, 3) = q(1:isize_field, 3)
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 3), txc(1, 4))
                txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                call Write_Visuals(plot_file, txc(:, 3:3))

                plot_file = 'PressureHydrostatic'//time_str(1:MaskSize)
                q = 0.0_wp
                call NSE_Pressure_Incompressible(q, s, txc(:, 2), txc(:, 3), txc(:, 6), txc(:, 7))
                call Write_Visuals(plot_file, txc(:, 2:2))

                plot_file = 'PressureHydrodynamic'//time_str(1:MaskSize)
                txc(1:isize_field, 1) = txc(1:isize_field, 1) - txc(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! ###################################################################
                ! Scalars
                ! ###################################################################
            case ('Scalars')
                do is = 1, inb_scal_array
                    write (str, *) is; plot_file = 'Scalar'//trim(adjustl(str))//time_str(1:MaskSize)

                    txc(1:isize_field, 1) = s(1:isize_field, is)
                    call Write_Visuals(plot_file, txc(:, 1:1))

                end do

                ! ###################################################################
                ! Scalars Derivatives
                ! ###################################################################
            case ('ScalarGradientVector')
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    plot_file = trim(adjustl(str))//'GradientVector'//time_str(1:MaskSize)
                    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s(1, is), txc(1, 1))
                    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s(1, is), txc(1, 2))
                    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s(1, is), txc(1, 3))
                    call Write_Visuals(plot_file, txc(:, 1:3))
                end do

            case ('ScalarGradient G_iG_i (Log)')
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    plot_file = 'Log'//trim(adjustl(str))//'Gradient'//time_str(1:MaskSize)
                    call FI_GRADIENT(imax, jmax, kmax, s(1, is), txc(1, 1), txc(1, 2))
                    txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                    call Write_Visuals(plot_file, txc(:, 1:1))

                end do

            case ('ScalarGradientEquation')
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    plot_file = trim(adjustl(str))//'Gradient'//time_str(1:MaskSize)
                    call FI_GRADIENT(imax, jmax, kmax, s(1, is), txc(1, 1), txc(1, 2))
                    call Write_Visuals(plot_file, txc(:, 1:1))

                    if (nse_diffusion == EQNS_NONE) then; diff = 0.0_wp
                    else; diff = visc/schmidt(is)
                    end if

                    plot_file = 'ScalarGradientProduction'//time_str(1:MaskSize)
                    call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s(1, is), q(1, 1), q(1, 2), q(1, 3), &
                                                txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                    call Write_Visuals(plot_file, txc(:, 1:1))

                    plot_file = trim(adjustl(str))//'GradientDiffusion'//time_str(1:MaskSize)
                    call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s(1, is), &
                                               txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                    txc(1:isize_field, 1) = diff*txc(1:isize_field, 1)
                    call Write_Visuals(plot_file, txc(:, 1:1))

                end do

                ! ###################################################################
                ! Velocity Derivatives; Vorticity
                ! ###################################################################
            case ('VorticityVector')
                call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call Write_Visuals(plot_file, txc(:, 1:3))

            case ('Enstrophy W_iW_i (Log)')
                plot_file = 'Log'//'Enstrophy'//time_str(1:MaskSize)
                call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'LogPotentialEnstrophy'//time_str(1:MaskSize)
                select case (nse_eqns)
                case (DNS_EQNS_BOUSSINESQ)
                    wrk1d(1:kmax, 1) = bbackground(1:kmax)
                    bbackground(1:kmax) = 0.0_wp
                    call Gravity_Source(gravityProps, imax, jmax, kmax, s, txc(:, 4))
                    bbackground(1:kmax) = wrk1d(1:kmax, 1)

                case (DNS_EQNS_ANELASTIC)
                    call Thermo_Anelastic_Buoyancy(imax, jmax, kmax, s, txc(:, 4))

                end select
                txc(1:isize_field, 4) = txc(1:isize_field, 4)/froude
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), txc(1, 4), txc(1, 1))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 4), txc(1, 2))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 4), txc(1, 3))
                call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
                txc(:, 1) = Tensor_Dot(txc(:, 1:3), txc(:, 4:6))
                txc(1:isize_field, 1) = log10(txc(1:isize_field, 1)*txc(1:isize_field, 1) + small_wp)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('EnstrophyEquation')
                plot_file = 'Enstrophy'//time_str(1:MaskSize)
                call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'EnstrophyProduction'//time_str(1:MaskSize)
                call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                             txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'EnstrophyDiffusion'//time_str(1:MaskSize)
                call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                            txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 1) = visc*txc(1:isize_field, 1)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! ###################################################################
                ! Velocity Derivatives; Strain
                ! ###################################################################
            case ('StrainTensor')
                call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call Write_Visuals(plot_file, txc(:, 1:6))

            case ('Strain 2S_ijS_ij (Log)')
                plot_file = 'LogStrain'//time_str(1:MaskSize)
                call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('StrainEquation')
                plot_file = 'Strain'//time_str(1:MaskSize)
                call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'StrainPressure'//time_str(1:MaskSize)
                if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                    call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
                else
                    txc(:, 1) = q(:, 6)     ! pressure
                end if
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 2) = 2.0_wp*txc(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 2:2))

                plot_file = 'StrainProduction'//time_str(1:MaskSize)
                call FI_STRAIN_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                          txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'StrainDiffusion'//time_str(1:MaskSize)
                call FI_STRAIN_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                         txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 1) = 2.0_wp*visc*txc(1:isize_field, 1)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! ###################################################################
                ! Velocity Derivatives; Invariants
                ! ###################################################################
            case ('VelocityGradientInvariants')
                plot_file = 'InvariantP'//time_str(1:MaskSize)
                call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'InvariantQ'//time_str(1:MaskSize)
                call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'InvariantR'//time_str(1:MaskSize)
                call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                    txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('HorizontalDivergence')
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), q(1, 1), txc(1, 1))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), q(1, 2), txc(1, 2))
                txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! ###################################################################
                ! Turbulence Quantities
                ! ###################################################################
            case ('Turbulent quantities')
                plot_file = 'Tke'//time_str(1:MaskSize)
                txc(1:isize_field, 1) = q(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 1))
                txc(1:isize_field, 2) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 2))
                txc(1:isize_field, 3) = q(1:isize_field, 3); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                txc(:, 4) = 0.5_wp*Tensor_Dot(txc(:, 1:3), txc(:, 1:3))
                call Write_Visuals(plot_file, txc(:, 4:4))

                plot_file = 'ReynoldsTensor'//time_str(1:MaskSize)
                txc(1:isize_field, 4) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
                txc(1:isize_field, 5) = txc(1:isize_field, 1)*txc(1:isize_field, 3)
                txc(1:isize_field, 6) = txc(1:isize_field, 2)*txc(1:isize_field, 3)
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 1)
                txc(1:isize_field, 2) = txc(1:isize_field, 2)*txc(1:isize_field, 2)
                txc(1:isize_field, 3) = txc(1:isize_field, 3)*txc(1:isize_field, 3)
                call Write_Visuals(plot_file, txc(:, 1:6))

                ! plot_file = 'LogDissipation'//time_str(1:MaskSize)
                ! call FI_DISSIPATION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                !                     txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
                ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*visc
                ! txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                ! call Write_Visuals(plot_file, txc(:,1:1))

                ! ###################################################################
                ! Buoyancy
                ! ###################################################################
            case ('Buoyancy')
                select case (nse_eqns)
                case (DNS_EQNS_BOUSSINESQ)
                    wrk1d(1:kmax, 1) = bbackground(1:kmax)
                    bbackground(1:kmax) = 0.0_wp
                    call Gravity_Source(gravityProps, imax, jmax, kmax, s, txc(1, 1))
                    bbackground(1:kmax) = wrk1d(1:kmax, 1)

                case (DNS_EQNS_ANELASTIC)
                    call Thermo_Anelastic_Buoyancy(imax, jmax, kmax, s, txc(:, 1))

                end select
                txc(1:isize_field, 1) = txc(1:isize_field, 1)/froude
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'Fwb'//time_str(1:MaskSize)     ! buoyancy flux along Oy
                txc(1:isize_field, 2) = txc(1:isize_field, 1)*q(1:isize_field, 3)
                call Write_Visuals(plot_file, txc(:, 2:2))

                plot_file = 'bPrime'//time_str(1:MaskSize)  ! buoyancy fluctuation
                call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(1, 1))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'Cwb'//time_str(1:MaskSize)     ! Covariance between b and w
                txc(1:isize_field, 2) = q(1:isize_field, 3); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 2))
                txc(1:isize_field, 2) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 2:2))

                plot_file = 'GradientRi'//time_str(1:MaskSize)
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(:, 1), txc(:, 2))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), q(:, 1), txc(:, 3))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), q(:, 2), txc(:, 4))
                txc(1:isize_field, 2) = abs(txc(1:isize_field, 2))/(txc(1:isize_field, 3)**2.0 + txc(1:isize_field, 4)**2.0 + small_wp)
                call Write_Visuals(plot_file, txc(:, 2:2))

                ! ###################################################################
                ! AirWater Anelastic
                ! ###################################################################
            case ('Anelastic')
                plot_file = 'RelativeHumidity'//time_str(1:MaskSize)
                call Thermo_Anelastic_RH(imax, jmax, kmax, s, txc(:, 1), wrk3d)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! plot_file = 'Dewpoint'//time_str(1:MaskSize)
                ! call Thermo_Anelastic_Weight_DewPoint(imax, jmax, kmax, s, txc(1, 1), wrk3d)
                ! call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'Theta'//time_str(1:MaskSize)
                call Thermo_Anelastic_Theta(imax, jmax, kmax, s, txc(:, 1), wrk3d)
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'ThetaV'//time_str(1:MaskSize)
                call Thermo_Anelastic_ThetaV(imax, jmax, kmax, s, txc(:, 1), wrk3d)
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'ThetaE'//time_str(1:MaskSize)
                call Thermo_Anelastic_ThetaE(imax, jmax, kmax, s, txc(:, 1), wrk3d)
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'ThetaL'//time_str(1:MaskSize)
                call Thermo_Anelastic_ThetaL(imax, jmax, kmax, s, txc(:, 1), wrk3d)
                call Write_Visuals(plot_file, txc(:, 1:1))

                ! ###################################################################
                ! Radiation
                ! ###################################################################
            case ('Infrared')
                do is = 1, inb_scal
                    if (infraredProps%active(is)) then
                        write (str, *) is; plot_file = 'Infrared'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Radiation_Infrared_Z(infraredProps, imax, jmax, kmax, fdm_Int0, s, &
                                                  txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4), txc(:, 5), txc(:, 6))
                        if (nse_eqns == DNS_EQNS_ANELASTIC) then
                            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, txc(:, 1))
                        end if
                        call Write_Visuals(plot_file, txc(:, 1:1))

                        write (str, *) is; plot_file = 'InfraredFlux'//trim(adjustl(str))//time_str(1:MaskSize)
                        txc(1:isize_field, 5) = txc(1:isize_field, 6) - txc(1:isize_field, 5)
                        call Write_Visuals(plot_file, txc(:, 1:5))

                        write (str, *) is; plot_file = 'InfraredFluxUp'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Write_Visuals(plot_file, txc(:, 1:6))
                    end if

                end do

                ! ###################################################################
                ! Microphysics
                ! ###################################################################
            case ('Sedimentation')
                do is = 1, inb_scal
                    if (infraredProps%active(is)) then
                        write (str, *) is; plot_file = 'Sedimentation'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Microphysics_Sedimentation(sedimentationProps, imax, jmax, kmax, is, g(3), s, &
                                                        txc(:, 1), txc(:, 2), txc(:, 3))
                        if (nse_eqns == DNS_EQNS_ANELASTIC) then
                            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, txc(:, 1))
                        end if
                        call Write_Visuals(plot_file, txc(:, 1:1))

                        write (str, *) is; plot_file = 'SedimentationFlux'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Write_Visuals(plot_file, txc(:, 1:3))

                    end if

                end do

            end select

        end do

    end do

    call TLab_Stop(0)

contains
    ! ###################################################################
    ! ###################################################################
    subroutine Visuals_Initialize()

        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        integer(wi) nz
#endif

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
        iv = 0
        iv = iv + 1; opt_name(iv) = 'VelocityX'
        iv = iv + 1; opt_name(iv) = 'VelocityY'
        iv = iv + 1; opt_name(iv) = 'VelocityZ'
        iv = iv + 1; opt_name(iv) = 'VelocityVector'
        iv = iv + 1; opt_name(iv) = 'VelocityMagnitude'
        iv = iv + 1; opt_name(iv) = 'Density'
        iv = iv + 1; opt_name(iv) = 'Temperature'
        iv = iv + 1; opt_name(iv) = 'Pressure'
        iv = iv + 1; opt_name(iv) = 'Scalars'
        iv = iv + 1; opt_name(iv) = 'ScalarGradientVector'
        iv = iv + 1; opt_name(iv) = 'ScalarGradient G_iG_i (Log)'
        iv = iv + 1; opt_name(iv) = 'ScalarGradientEquation'
        iv = iv + 1; opt_name(iv) = 'VorticityVector'
        iv = iv + 1; opt_name(iv) = 'Enstrophy W_iW_i (Log)'
        iv = iv + 1; opt_name(iv) = 'EnstrophyEquation'
        iv = iv + 1; opt_name(iv) = 'StrainTensor'
        iv = iv + 1; opt_name(iv) = 'Strain 2S_ijS_ij (Log)'
        iv = iv + 1; opt_name(iv) = 'StrainEquation'
        iv = iv + 1; opt_name(iv) = 'VelocityGradientInvariants'
        iv = iv + 1; opt_name(iv) = 'HorizontalDivergence'
        iv = iv + 1; opt_name(iv) = 'Turbulent quantities'
        iv = iv + 1; opt_name(iv) = 'Buoyancy'
        iv = iv + 1; opt_name(iv) = 'Sedimentation'
        iv = iv + 1; opt_name(iv) = 'Infrared'
        iv = iv + 1; opt_name(iv) = 'Anelastic'
        if (iv > iopt_size_max) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Increase number of options.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! 17, '. Relative humidity'
        ! 19, '. Analysis of B and V'
        ! 20, '. Total Stress Tensor'

        call ScanFile_Char(bakfile, ifile, block, 'ParamVisuals', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, '(A)') 'Option?'
            do is = 1, iv
                write (*, '(I2,A)') is, '. '//trim(adjustl(opt_name(is)))
            end do
            read (*, '(A512)') sRes
#endif
        end if
        iopt_size = iopt_size_max
        call LIST_INTEGER(sRes, iopt_size, opt_vec)

        if (opt_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing ParamVisuals input.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! -------------------------------------------------------------------
        call ScanFile_Char(bakfile, ifile, block, 'Subdomain', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Subdomain limits (full domain, default)?'
            read (*, '(A64)') sRes
#endif
        end if
        is = 6
        call LIST_INTEGER(sRes, is, subdomain)

        if (is < 6) then ! default
            subdomain(1) = 1; subdomain(2) = x%size
            subdomain(3) = 1; subdomain(4) = y%size
            subdomain(5) = 1; subdomain(6) = z%size
        end if

        ! -------------------------------------------------------------------
        call ScanFile_Int(bakfile, ifile, block, 'Format', '-1', opt_format)
        if (opt_format == -1) then
#ifdef USE_MPI
#else
            write (*, *) 'File Format ?'
            write (*, *) ' 1. Raw, single precision, no header (default)'
            write (*, *) ' 2. General restart format'
            read (*, '(A64)') sRes
#endif
        end if
        if (len_trim(adjustl(sRes)) > 0) then
            if (trim(adjustl(sRes)) == 'single') then; opt_format = FORMAT_SINGLE
            else if (trim(adjustl(sRes)) == 'general') then; opt_format = FORMAT_GENERAL
            else
                read (sRes, *) opt_format
            end if
        end if

        if (opt_format < 0) opt_format = FORMAT_SINGLE ! default is single precission, no header

!     ! -------------------------------------------------------------------
!     ! Defining gate levels for conditioning
!     ! -------------------------------------------------------------------
!     opt_cond = 0 ! default values
!     opt_cond_relative = 0
!     igate_size = 0

!     do iv = 1, iopt_size
!         if (opt_vec(iv) == iscal_offset + 11) then
! #include "dns_read_partition.h"
!             if (opt_cond > 1) inb_txc = max(inb_txc, 5)
!             exit
!         end if
!     end do

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        iread_flow = .false.
        iread_scal = .false.
        inb_txc = 0

        do iv = 1, iopt_size
            select case (trim(adjustl(opt_name(opt_vec(iv)))))
            case ('VelocityX')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('VelocityY')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('VelocityZ')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('VelocityVector')
                iread_flow = .true.; inb_txc = max(inb_txc, 3)
            case ('VelocityMagnitude')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('Density')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 1)
            case ('Temperature')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 1)
            case ('Pressure')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7)
            case ('Scalars')
                iread_scal = .true.; inb_txc = max(inb_txc, 1)
            case ('ScalarGradientVector')
                iread_scal = .true.; inb_txc = max(inb_txc, 3)
            case ('ScalarGradient G_iG_i (Log)')
                iread_scal = .true.; inb_txc = max(inb_txc, 2)
            case ('ScalarGradientEquation')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 6)
            case ('VorticityVector')
                iread_flow = .true.; inb_txc = max(inb_txc, 4)
            case ('Enstrophy W_iW_i (Log)')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7)
            case ('EnstrophyEquation')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('StrainTensor')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('Strain 2S_ijS_ij (Log)')
                iread_flow = .true.; inb_txc = max(inb_txc, 3)
            case ('StrainEquation')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('VelocityGradientInvariants')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('HorizontalDivergence')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 2)
            case ('Turbulent quantities')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('Buoyancy')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 4)
            case ('Anelastic')
                iread_scal = .true.; inb_txc = max(inb_txc, 1)
            case ('Infrared')
                iread_scal = .true.; inb_txc = max(inb_txc, 6)
            case ('Sedimentation')
                iread_scal = .true.; inb_txc = max(inb_txc, 3)
            end select
        end do

        ! #######################################################################
#ifdef USE_MPI

        io_subarrays(:)%active = .false.
        io_subarrays(:)%offset = 0
        io_subarrays(:)%precision = IO_TYPE_SINGLE

        nz = subdomain(6) - subdomain(5) + 1

        ! ###################################################################
        ! Saving full vertical xOz planes; using subdomain(3) to define the plane
        if (ims_pro_j == ((subdomain(3) - 1)/jmax)) io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%communicator = ims_comm_x
        io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%subarray = IO_Create_Subarray_XOZ(imax, nz, MPI_REAL4)

        ! Saving full vertical yOz planes; using subiddomain(1) to define the plane
        if (ims_pro_i == ((subdomain(1) - 1)/imax)) io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%communicator = ims_comm_y
        io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%subarray = IO_Create_Subarray_YOZ(jmax, nz, MPI_REAL4)

        ! Saving full blocks xOy planes
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%communicator = MPI_COMM_WORLD
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%subarray = IO_Create_Subarray_XOY(imax, jmax, nz, MPI_REAL4)

#else
        io_subarrays(:)%offset = 0
        io_subarrays(:)%precision = IO_TYPE_SINGLE

#endif

        return
    end subroutine Visuals_Initialize

!########################################################################
!########################################################################
    subroutine Write_Visuals(fname, field)
        use Reductions, only: Reduce_Block_InPlace

        character(len=*) fname
        real(wp), intent(inout) :: field(:, :)

        ! -------------------------------------------------------------------
        integer(wi) nx, ny, nz, ifield, i
        integer(wi) sizes(5)
        character*32 varname(16)
        integer subarray_plan
        integer nfield

        ! ###################################################################
        ! I think this first block should go into Visuals_Initialize
        nfield = size(field, 2)

        nx = subdomain(2) - subdomain(1) + 1
        ny = subdomain(4) - subdomain(3) + 1
        nz = subdomain(6) - subdomain(5) + 1

        sizes(5) = nfield
        sizes(1) = size(field, 1)           ! array size
        sizes(2) = 1                        ! lower bound
        if (subdomain(2) - subdomain(1) + 1 == x%size .and. &
            subdomain(4) - subdomain(3) + 1 == 1) then              ! xOz plane
            subarray_plan = IO_SUBARRAY_VISUALS_XOZ
            sizes(3) = imax*nz              ! upper bound
            sizes(4) = 1                    ! stride
            nx = imax

        else if (subdomain(4) - subdomain(3) + 1 == y%size .and. &
                 subdomain(2) - subdomain(1) + 1 == 1) then         ! yOz plane
            subarray_plan = IO_SUBARRAY_VISUALS_YOZ
            sizes(3) = jmax*nz              ! upper bound
            sizes(4) = 1                    ! stride
            ny = jmax

        else if (subdomain(2) - subdomain(1) + 1 == x%size .and. &
                 subdomain(4) - subdomain(3) + 1 == y%size) then    ! xOy blocks
            subarray_plan = IO_SUBARRAY_VISUALS_XOY
            sizes(3) = imax*jmax*nz         ! upper bound
            sizes(4) = 1                    ! stride
            nx = imax
            ny = jmax

        else
#ifdef USE_MPI
            call TLab_Write_ASCII(efile, __FILE__//'. Invalid subdomain in parallel mode.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
#else
            sizes(3) = nx*ny*nz             ! upper bound
            sizes(4) = 1                    ! stride
#endif

        end if

        ! ###################################################################
        select case (opt_format)
        case (FORMAT_GENERAL)
            if (nfield > 1) then ! IO_Write_Fields expects field to be aligned by isize_field (instead of isize_txc_field)
                do ifield = 2, nfield
                    do i = 1, isize_field
                        field((ifield - 1)*isize_field + i, 1) = field(i, ifield)
                    end do
                end do
            end if
            call IO_Write_Fields(fname, imax, jmax, kmax, itime, nfield, field)

        case (FORMAT_SINGLE)
            do ifield = 1, nfield
                call Reduce_Block_InPlace(imax, jmax, kmax, subdomain(1), subdomain(3), subdomain(5), nx, ny, nz, field(:, ifield))
            end do

            varname = ''
            if (nfield > 1) then
                do ifield = 1, nfield; write (varname(ifield), *) ifield; varname(ifield) = trim(adjustl(varname(ifield)))
                end do
            end if
            call IO_Write_Subarray(io_subarrays(subarray_plan), fname, varname, field, sizes)

        end select

        return
    end subroutine Write_Visuals

end program VISUALS
