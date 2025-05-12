#include "tlab_error.h"

program VISUALS
    use TLab_Constants, only: wp, wi, small_wp, MAX_PARS, fmt_r
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_txc, inb_flow, inb_scal, isize_field, isize_txc_field
    use TLab_Arrays
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro, ims_pro_i, ims_pro_j, ims_comm_x, ims_comm_y
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: g, FDM_Initialize
    ! use FDM, only: fdm_Int0
    ! use Thermodynamics, only: NSP, THERMO_SPNAME
    ! use Thermodynamics, only: Thermodynamics_Initialize_Parameters
    ! use Thermodynamics, only: imixture, MIXT_TYPE_NONE, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR
    use NavierStokes
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize, gravityProps, Gravity_Source
    ! use Rotation, only: Rotation_Initialize
    ! use Thermo_Anelastic
    ! use THERMO_AIRWATER
    ! use Radiation
    ! use Microphysics
    ! use Chemistry
    ! use LargeScaleForcing, only: LargeScaleForcing_Initialize
    ! use PARTICLE_VARS
    ! use PARTICLE_ARRAYS
    ! use PARTICLE_PROCS
    use OPR_Partial
    use OPR_Fourier
    use OPR_Elliptic
    use OPR_Burgers, only: OPR_Burgers_Initialize
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN

    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! Iteration # to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 20    ! Maximum # of options
    integer(wi) iopt_size, iv
    integer(wi) opt_vec(iopt_size_max)
    ! real(wp) opt_vec2(iopt_size_max)

    integer :: opt_format                           ! File format
    integer, parameter :: FORMAT_SINGLE = 1         ! Single precision, no headers
    integer, parameter :: FORMAT_GENERAL = 2        ! General IO format

    integer(wi) subdomain(6)                        ! Subdomain to be saved

    character(len=32) bakfile, block
    character(len=128) eStr
    character(len=512) sRes

    ! character*32 fname
    character*32 flow_file, scal_file, plot_file
    character*64 str

    character(len=32) time_str
    integer, parameter :: MaskSize = 6

    integer, parameter :: IO_SUBARRAY_VISUALS_XOY = 1
    integer, parameter :: IO_SUBARRAY_VISUALS_YOZ = 2
    integer, parameter :: IO_SUBARRAY_VISUALS_XOZ = 3
    type(io_subarray_dt) :: io_subarrays(3)

    ! integer(wi) opt_cond, opt_cond_scal, opt_cond_relative
    integer(wi) ij, is!, ig
    integer(wi) iscal_offset, idummy
    logical iread_flow, iread_scal !, iread_part
    real(wp) diff
    real(wp) params(MAX_PARS)

    ! ! Gates for the definition of the intermittency function (partition of the fields)
    ! integer(wi), parameter :: igate_size_max = 8
    ! integer(wi) igate_size
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
    ! call Particle_Initialize_Parameters(ifile)

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    ! call Thermodynamics_Initialize_Parameters(ifile)

    call Gravity_Initialize(ifile)
    ! call Rotation_Initialize(ifile)
    ! call Radiation_Initialize(ifile)
    ! call Microphysics_Initialize(ifile)
    ! call LargeScaleForcing_Initialize(ifile)
    ! call Chemistry_Initialize(ifile)
    ! call Particle_Initialize_Parameters(ifile)

    call TLab_Consistency_Check()

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
    iscal_offset = 9    ! define iscal_offset
    ! if (imixture == MIXT_TYPE_NONE) then
    !     iscal_offset = 9 + NSP
    ! end if

    call ScanFile_Char(bakfile, ifile, block, 'ParamVisuals', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_INTEGER(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, '(A)') 'Option?'
        ! write (*, '(A)') ' 0. Grid'
        write (*, '(A)') ' 1. VelocityX'
        write (*, '(A)') ' 2. VelocityY'
        write (*, '(A)') ' 3. VelocityZ'
        write (*, '(A)') ' 4. VelocityVector'
        write (*, '(A)') ' 5. Velocity V_iV_i'
        write (*, '(A)') ' 6. Density'
        write (*, '(A)') ' 7. Temperature'
        write (*, '(A)') ' 8. Pressure'
        write (*, '(A)') ' 9. Scalars'
        ! if (imixture /= MIXT_TYPE_NONE) then
        !     do is = 1, NSP
        !         write (*, '(I2,A)') 9 + is, '. '//THERMO_SPNAME(is)
        !     end do
        ! end if

        write (*, '(I2,A)') iscal_offset + 1, '. ScalarGradientVector'
        write (*, '(I2,A)') iscal_offset + 2, '. ScalarGradient G_iG_i (Log)'
        write (*, '(I2,A)') iscal_offset + 3, '. ScalarGradientEquation'
        write (*, '(I2,A)') iscal_offset + 4, '. VorticityVector'
        write (*, '(I2,A)') iscal_offset + 5, '. Enstrophy W_iW_i (Log)'
        write (*, '(I2,A)') iscal_offset + 6, '. EnstrophyEquation'
        write (*, '(I2,A)') iscal_offset + 7, '. StrainTensor'
        write (*, '(I2,A)') iscal_offset + 8, '. Strain 2S_ijS_ij (Log)'
        write (*, '(I2,A)') iscal_offset + 9, '. StrainEquation'
        write (*, '(I2,A)') iscal_offset + 10, '. VelocityGradientInvariants'
        write (*, '(I2,A)') iscal_offset + 11, '. Space partition'
        write (*, '(I2,A)') iscal_offset + 12, '. Buoyancy'
        write (*, '(I2,A)') iscal_offset + 14, '. HorizontalDivergence'
        write (*, '(I2,A)') iscal_offset + 15, '. Turbulent quantities'
        write (*, '(I2,A)') iscal_offset + 16, '. Radiative forcing'
        write (*, '(I2,A)') iscal_offset + 17, '. Relative humidity'
        write (*, '(I2,A)') iscal_offset + 18, '. Particle Density'
        write (*, '(I2,A)') iscal_offset + 19, '. Analysis of B and V'
        write (*, '(I2,A)') iscal_offset + 20, '. Total Stress Tensor'
        read (*, '(A512)') sRes
#endif
    end if
    iopt_size = iopt_size_max
    call LIST_INTEGER(sRes, iopt_size, opt_vec)

    if (opt_vec(1) < 0) then ! Check
        call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing input [PostProcessing.ParamVisuals].')
        call TLab_Stop(DNS_ERROR_INVALOPT)
    end if

    ! ###################################################################
    ! Initialization of array sizes
    ! ###################################################################
    iread_flow = .false.
    iread_scal = .false.
    inb_txc = 0

    do iv = 1, iopt_size
        if (opt_vec(iv) == 1) then; iread_flow = .true.; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 2) then; iread_flow = .true.; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 3) then; iread_flow = .true.; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 4) then; iread_flow = .true.; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == 5) then; iread_flow = .true.; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) == 6) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == 7) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == 8) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7); end if
        if (opt_vec(iv) == 8) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 10); end if
        if (opt_vec(iv) == 9) then; iread_scal = .true.; inb_txc = max(inb_txc, 1); end if
        if (opt_vec(iv) > 9 .and. opt_vec(iv) <= iscal_offset) then
            iread_scal = .true.; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 1) then; iread_scal = .true.; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == iscal_offset + 2) then; iread_scal = .true.; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == iscal_offset + 3) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 4) then; iread_flow = .true.; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 5) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7); end if
        if (opt_vec(iv) == iscal_offset + 6) then; iread_flow = .true.; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 7) then; iread_flow = .true.; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 8) then; iread_flow = .true.; inb_txc = max(inb_txc, 3); end if
        if (opt_vec(iv) == iscal_offset + 9) then; iread_flow = .true.; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 10) then; iread_flow = .true.; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 12) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 14) then; iread_flow = .true.; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 15) then; iread_flow = .true.; inb_txc = max(inb_txc, 6); end if
        if (opt_vec(iv) == iscal_offset + 16) then; iread_scal = .true.; inb_txc = max(inb_txc, 4); end if
        if (opt_vec(iv) == iscal_offset + 17) then; iread_scal = .true.; inb_txc = max(inb_txc, 2); end if
        ! if (opt_vec(iv) == iscal_offset + 18) then; iread_part = .true.; inb_txc = max(inb_txc, 2); end if
        if (opt_vec(iv) == iscal_offset + 19) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7); end if
        if (opt_vec(iv) == iscal_offset + 20) then; iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 10); end if
    end do

    ! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, ifile, block, 'Subdomain', '-1', sRes)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'Subdomain limits (full domain, default)?'
        read (*, '(A64)') sRes
#endif
    end if
    idummy = 6
    call LIST_INTEGER(sRes, idummy, subdomain)

    if (idummy < 6) then ! default
        subdomain(1) = 1; subdomain(2) = x%size
        subdomain(3) = 1; subdomain(4) = y%size
        subdomain(5) = 1; subdomain(6) = z%size
    end if

    ! -------------------------------------------------------------------
    opt_format = FORMAT_SINGLE ! default values

    call ScanFile_Char(bakfile, ifile, block, 'Format', '-1', sRes)

    if (sRes == '-1') then
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

    ! allocate (gate(isize_field))

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call OPR_Partial_Initialize()

    call TLab_Initialize_Background(ifile)

    call OPR_Burgers_Initialize(ifile)

    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

#ifdef USE_MPI
    call Create_IO_Subarrays()
#else
    io_subarrays(:)%offset = 0
    io_subarrays(:)%precision = IO_TYPE_SINGLE
#endif

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))//'.'
        call TLab_Write_ASCII(lfile, sRes)

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

        ! call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

        write (sRes, fmt_r) rtime; sRes = 'Physical time '//trim(adjustl(sRes))
        call TLab_Write_ASCII(lfile, sRes)

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

            ! ###################################################################
            ! Velocities
            ! ###################################################################
            if (opt_vec(iv) == 1) then
                txc(1:isize_field, 1) = q(1:isize_field, 1)
                plot_file = 'VelocityX'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            else if (opt_vec(iv) == 2) then
                txc(1:isize_field, 1) = q(1:isize_field, 2)
                plot_file = 'VelocityY'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            else if (opt_vec(iv) == 3) then
                txc(1:isize_field, 1) = q(1:isize_field, 3)
                plot_file = 'VelocityZ'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            else if (opt_vec(iv) == 4) then
                txc(1:isize_field, 1:3) = q(1:isize_field, 1:3)
                plot_file = 'VelocityVector'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 3, subdomain, txc(1, 1))

            else if (opt_vec(iv) == 5) then
                txc(1:isize_field, 1) = sqrt(q(1:isize_field, 1)**2 + q(1:isize_field, 2)**2 + q(1:isize_field, 3)**2)
                plot_file = 'VelocityMagnitude'//time_str(1:MaskSize)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            end if

            ! ###################################################################
            ! Thermodynamic state
            ! ###################################################################
            if (opt_vec(iv) == 6) then
                plot_file = 'Density'//time_str(1:MaskSize)
                select case (nse_eqns)
                case (DNS_EQNS_COMPRESSIBLE)
                    txc(1:isize_field, 1) = q(1:isize_field, 5)

                case (DNS_EQNS_BOUSSINESQ)
                    ! wrk1d(1:jmax, 1) = 0.0_wp
                    ! call Gravity_Source(gravity, imax, jmax, kmax, s, txc(1, 1), wrk1d)
                    ! dummy = 1.0_wp/froude
                    ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy + 1.0_wp

                case (DNS_EQNS_ANELASTIC)
                    ! call Thermo_Anelastic_DENSITY(imax, jmax, kmax, s, txc(1, 1), wrk3d)

                end select
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            else if (opt_vec(iv) == 7) then
                plot_file = 'Temperature'//time_str(1:MaskSize)
                select case (nse_eqns)
                case (DNS_EQNS_COMPRESSIBLE)
                    txc(1:isize_field, 1) = q(1:isize_field, 7)

                case (DNS_EQNS_BOUSSINESQ)
                case (DNS_EQNS_ANELASTIC)
                    ! call Thermo_Anelastic_TEMPERATURE(imax, jmax, kmax, s, txc(1, 1))

                end select
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            else if (opt_vec(iv) == 8) then
                plot_file = 'Pressure'//time_str(1:MaskSize)
                select case (nse_eqns)
                case (DNS_EQNS_COMPRESSIBLE)
                    txc(1:isize_field, 1) = q(1:isize_field, 6)
                case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
                    ! call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                end select
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                ! Additional pressure stuff...
                plot_file = 'PressureGradientPower'//time_str(1:MaskSize)
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), txc(1, 1), txc(1, 2))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 1), txc(1, 3))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 1), txc(1, 4))
                txc(1:isize_field, 2) = -(txc(1:isize_field, 2)*q(1:isize_field, 1) &
                                          + txc(1:isize_field, 3)*q(1:isize_field, 2) &
                                          + txc(1:isize_field, 4)*q(1:isize_field, 3))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

                call TLab_Write_ASCII(lfile, 'Computing pressure-strain correlation...')
                txc(1:isize_field, 2) = txc(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(1, 2))

                plot_file = 'PressureStrainX'//time_str(1:MaskSize)
                txc(1:isize_field, 3) = q(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), txc(1, 3), txc(1, 4))
                txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 3))

                plot_file = 'PressureStrainY'//time_str(1:MaskSize)
                txc(1:isize_field, 3) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 3), txc(1, 4))
                txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 3))

                plot_file = 'PressureStrainZ'//time_str(1:MaskSize)
                txc(1:isize_field, 3) = q(1:isize_field, 3); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 3), txc(1, 4))
                txc(1:isize_field, 3) = txc(1:isize_field, 2)*txc(1:isize_field, 4)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 3))

                plot_file = 'PressureHydrostatic'//time_str(1:MaskSize)
                q = 0.0_wp
                ! call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

                plot_file = 'PressureHydrodynamic'//time_str(1:MaskSize)
                txc(1:isize_field, 1) = txc(1:isize_field, 1) - txc(1:isize_field, 2)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            end if

            ! ###################################################################
            ! Scalars
            ! ###################################################################
            if (opt_vec(iv) == 9) then ! All prognostic scalars
                do is = 1, inb_scal_array
                    write (str, *) is; plot_file = 'Scalar'//trim(adjustl(str))//time_str(1:MaskSize)
                    txc(1:isize_field, 1) = s(1:isize_field, is)
                    call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
                end do

            else if (opt_vec(iv) > 9 .and. opt_vec(iv) <= iscal_offset) then ! Individual and diagnostic scalars
                ! if (imixture == MIXT_TYPE_AIRWATER) then ! s(1,inb_scal+1) contains liquid mass fraction
                !     if (opt_vec(iv) == 10) then ! vapor water mass fraction
                !         plot_file = trim(adjustl(THERMO_SPNAME(1)))//time_str(1:MaskSize)
                !         txc(1:isize_field, 1) = s(1:isize_field, inb_scal) - s(1:isize_field, inb_scal + 1)
                !         call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                !     else if (opt_vec(iv) == 11) then ! air mass fraction
                !         plot_file = trim(adjustl(THERMO_SPNAME(2)))//time_str(1:MaskSize)
                !         txc(1:isize_field, 1) = 1.0_wp - s(1:isize_field, inb_scal)
                !         call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                !     else if (opt_vec(iv) == 12) then ! liquid mass fraction
                !         plot_file = trim(adjustl(THERMO_SPNAME(3)))//time_str(1:MaskSize)
                !         call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, s(1, inb_scal + 1))

                !     end if

                ! else ! Plot the chosen species
                !     is = opt_vec(iv) - 9; plot_file = trim(adjustl(THERMO_SPNAME(is)))//time_str(1:MaskSize)
                !     txc(1:isize_field, 1) = s(1:isize_field, is)
                !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                ! end if

            end if

            ! ###################################################################
            ! Scalar Derivatives
            ! ###################################################################
            if (opt_vec(iv) >= iscal_offset + 1 .and. opt_vec(iv) <= iscal_offset + 3) then
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    if (opt_vec(iv) == iscal_offset + 1) then
                        plot_file = trim(adjustl(str))//'GradientVector'//time_str(1:MaskSize)
                        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s(1, is), txc(1, 1))
                        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s(1, is), txc(1, 2))
                        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s(1, is), txc(1, 3))
                        call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 3, subdomain, txc(1, 1))
                    end if

                    if (opt_vec(iv) == iscal_offset + 2 .or. opt_vec(iv) == iscal_offset + 3) then ! Gradient magnitude
                        plot_file = trim(adjustl(str))//'Gradient'//time_str(1:MaskSize)
                        call FI_GRADIENT(imax, jmax, kmax, s(1, is), txc(1, 1), txc(1, 2))
                        if (opt_vec(iv) == iscal_offset + 2) then
                            plot_file = 'Log'//trim(adjustl(plot_file))
                            txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                        end if
                        call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
                    end if

                    if (opt_vec(iv) == iscal_offset + 3 .and. is <= inb_scal) then ! Scalar gradient equation
                        if (nse_diffusion == EQNS_NONE) then; diff = 0.0_wp
                        else; diff = visc/schmidt(is)
                        end if

                        call TLab_Write_ASCII(lfile, 'Computing scalar gradient production...')
                        plot_file = 'ScalarGradientProduction'//time_str(1:MaskSize)
                        call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s(1, is), q(1, 1), q(1, 2), q(1, 3), &
                                                    txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                        call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                        call TLab_Write_ASCII(lfile, 'Computing scalar gradient diffusion...')
                        plot_file = trim(adjustl(str))//'GradientDiffusion'//time_str(1:MaskSize)
                        call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s(1, is), &
                                                   txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                        txc(1:isize_field, 1) = diff*txc(1:isize_field, 1)
                        call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                    end if

                end do

            end if

            ! ###################################################################
            ! Velocity Derivatives
            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 4) then ! VorticityVector
                plot_file = 'VorticityVector'//time_str(1:MaskSize)
                call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 3, subdomain, txc(1, 1))
            end if

            if (opt_vec(iv) == iscal_offset + 5 .or. opt_vec(iv) == iscal_offset + 6) then ! Enstrophy
                plot_file = 'Enstrophy'//time_str(1:MaskSize)
                call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                if (opt_vec(iv) == iscal_offset + 5) then ! Natural log
                    plot_file = 'Log'//trim(adjustl(plot_file))
                    txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                end if
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                ! plot_file = 'LogPotentialEnstrophy'//time_str(1:MaskSize)
                ! if (buoyancy%type == EQNS_BOD_EXPLICIT) then
                !     call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, txc(1, 4))
                ! else
                !     wrk1d(1:jmax, 1) = 0.0_wp
                !     call Gravity_Source(buoyancy, imax, jmax, kmax, s, txc(1, 4), wrk1d)
                ! end if
                ! dummy = 1.0_wp/froude
                ! txc(1:isize_field, 4) = txc(1:isize_field, 4)*dummy
                ! call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), txc(1, 4), txc(1, 1))
                ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 4), txc(1, 2))
                ! call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), txc(1, 4), txc(1, 3))
                ! call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
                ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4) &
                !                         + txc(1:isize_field, 2)*txc(1:isize_field, 5) &
                !                         + txc(1:isize_field, 3)*txc(1:isize_field, 6)
                ! txc(1:isize_field, 1) = log10(txc(1:isize_field, 1)*txc(1:isize_field, 1) + small_wp)
                ! call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
            end if

            if (opt_vec(iv) == iscal_offset + 6) then ! EnstrophyEquation
                call TLab_Write_ASCII(lfile, 'Computing enstrophy production...')
                plot_file = 'EnstrophyProduction'//time_str(1:MaskSize)
                call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                             txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                call TLab_Write_ASCII(lfile, 'Computing enstrophy diffusion...')
                plot_file = 'EnstrophyDiffusion'//time_str(1:MaskSize)
                call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                            txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 1) = visc*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            end if

            ! -------------------------------------------------------------------
            if (opt_vec(iv) == iscal_offset + 7) then ! Strain Tensor
                plot_file = 'StrainTensor'//time_str(1:MaskSize)
                call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 6, subdomain, txc(1, 1))
            end if

            if (opt_vec(iv) == iscal_offset + 8 .or. opt_vec(iv) == iscal_offset + 9) then ! Strain
                plot_file = 'Strain'//time_str(1:MaskSize)
                call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                if (opt_vec(iv) == iscal_offset + 8) then ! Natural log
                    plot_file = 'Log'//trim(adjustl(plot_file))
                    txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                end if
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
            end if

            if (opt_vec(iv) == iscal_offset + 9) then ! StrainEquation (I need the pressure)
                call TLab_Write_ASCII(lfile, 'Computing strain pressure...')
                plot_file = 'StrainPressure'//time_str(1:MaskSize)
                if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                    call TLab_Write_ASCII(efile, 'VISUALS. Strain eqn for incompressible undeveloped.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)
                else
                    txc(:, 6) = q(:, 6)
                end if
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 6), &
                                        txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                call TLab_Write_ASCII(lfile, 'Computing strain production...')
                plot_file = 'StrainProduction'//time_str(1:MaskSize)
                call FI_STRAIN_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                          txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                call TLab_Write_ASCII(lfile, 'Computing strain diffusion...')
                plot_file = 'StrainDiffusion'//time_str(1:MaskSize)
                call FI_STRAIN_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                         txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                txc(1:isize_field, 1) = 2.0_wp*visc*txc(1:isize_field, 1)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            end if

            ! -------------------------------------------------------------------
            if (opt_vec(iv) == iscal_offset + 10) then ! Invariants
                plot_file = 'InvariantP'//time_str(1:MaskSize)
                call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                plot_file = 'InvariantQ'//time_str(1:MaskSize)
                call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                plot_file = 'InvariantR'//time_str(1:MaskSize)
                call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                    txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            end if

            ! ###################################################################
            ! Partition
            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 11) then
                call TLab_Write_ASCII(efile, 'VISUALS. Partition undevelop.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            ! ! ###################################################################
            ! ! Buoyancy
            ! ! ###################################################################
            ! if (opt_vec(iv) == iscal_offset + 12) then
            !     plot_file = 'Buoyancy'//time_str(1:MaskSize)
            !     if (gravityProps%type == EQNS_BOD_EXPLICIT) then
            !         call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, txc(1, 1))
            !     else
            !         wrk1d(1:jmax, 1) = 0.0_wp
            !         call Gravity_Source(gravityProps, imax, jmax, kmax, s, txc(1, 1), wrk1d)
            !     end if
            !     dummy = 1.0_wp/froude
            !     txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            !     plot_file = 'Fvb'//time_str(1:MaskSize)     ! buoyancy flux along Oy
            !     txc(1:isize_field, 2) = txc(1:isize_field, 1)*q(1:isize_field, 2)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

            !     plot_file = 'bPrime'//time_str(1:MaskSize)  ! buoyancy fluctuation
            !     call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(1, 1))
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            !     plot_file = 'Cvb'//time_str(1:MaskSize)     ! Covariance between b and v
            !     txc(1:isize_field, 2) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 2))
            !     txc(1:isize_field, 2) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

            !     plot_file = 'LogBuoyancySource'//time_str(1:MaskSize)
            !     if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            !         call THERMO_AIRWATER_LINEAR_SOURCE(imax*jmax*kmax, s, txc(1, 1), txc(1, 2), txc(1, 3))
            !         call FI_GRADIENT(imax, jmax, kmax, txc(1, 1), txc(1, 2), txc(1, 4))
            !         dummy = gravityProps%parameters(inb_scal_array)
            !         txc(1:isize_field, 2) = txc(1:isize_field, 2)*txc(1:isize_field, 3)*dummy
            !     else
            !         call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 1), txc(1, 2))
            !         call Gravity_Source_Source(gravityProps, imax, jmax, kmax, s, txc(1, 1), txc(1, 2))
            !     end if
            !     dummy = visc/schmidt(1)/froude
            !     txc(1:isize_field, 1) = txc(1:isize_field, 2)*dummy
            !     txc(1:isize_field, 1) = log10(abs(txc(1:isize_field, 1)) + small_wp)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            ! end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 14) then
                plot_file = 'HorizontalDivergence'//time_str(1:MaskSize)
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), q(1, 1), txc(1, 2))
                call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), q(1, 3), txc(1, 1))
                txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
            end if

            ! ###################################################################
            if (opt_vec(iv) == iscal_offset + 15) then ! Turbulent quantities
                ! plot_file = 'LogDissipation'//time_str(1:MaskSize)
                ! call FI_DISSIPATION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                !                     txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
                ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*visc
                ! txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                ! call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

                plot_file = 'Tke'//time_str(1:MaskSize)
                txc(1:isize_field, 1) = q(1:isize_field, 1); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 1))
                txc(1:isize_field, 2) = q(1:isize_field, 2); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 2))
                txc(1:isize_field, 3) = q(1:isize_field, 3); call FI_FLUCTUATION_INPLACE(imax, jmax, kmax, txc(:, 3))
                txc(1:isize_field, 4) = 0.5_wp*(txc(1:isize_field, 1)**2 + txc(1:isize_field, 2)**2 + txc(1:isize_field, 3)**2)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 4))

                plot_file = 'ReynoldsTensor'//time_str(1:MaskSize)
                txc(1:isize_field, 4) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
                txc(1:isize_field, 5) = txc(1:isize_field, 1)*txc(1:isize_field, 3)
                txc(1:isize_field, 6) = txc(1:isize_field, 2)*txc(1:isize_field, 3)
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 1)
                txc(1:isize_field, 2) = txc(1:isize_field, 2)*txc(1:isize_field, 2)
                txc(1:isize_field, 3) = txc(1:isize_field, 3)*txc(1:isize_field, 3)
                call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 6, subdomain, txc(1, 1))

            end if

            ! ! ###################################################################
            ! if (opt_vec(iv) == iscal_offset + 16) then
            !     do is = 1, inb_scal

            !         if (infraredProps%active(is)) then
            !             write (str, *) is; plot_file = 'Radiation'//trim(adjustl(str))//time_str(1:MaskSize)
            !             call Radiation_Infrared_Y(infraredProps, imax, jmax, kmax, fdm_Int0, s, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
            !             call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
            !         end if

            !     end do
            ! end if

            ! ! ###################################################################
            ! if (opt_vec(iv) == iscal_offset + 17) then
            !     plot_file = 'RelativeHumidity'//time_str(1:MaskSize)
            !     call Thermo_Anelastic_RELATIVEHUMIDITY(imax, jmax, kmax, s, txc(1, 1), wrk3d)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))
            ! end if

            ! ! ###################################################################
            ! if (opt_vec(iv) == iscal_offset + 19) then
            !     plot_file = 'LaplacianV'//time_str(1:MaskSize)
            !     call OPR_Partial_Z(OPR_P2, imax, jmax, kmax, g(3), q(1, 2), txc(1, 4), txc(1, 5))
            !     call OPR_Partial_Y(OPR_P2, imax, jmax, kmax, g(2), q(1, 2), txc(1, 3), txc(1, 5))
            !     call OPR_Partial_X(OPR_P2, imax, jmax, kmax, g(1), q(1, 2), txc(1, 2), txc(1, 5))
            !     txc(1:isize_field, 2) = txc(1:isize_field, 2) + txc(1:isize_field, 3) + txc(1:isize_field, 4)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

            !     plot_file = 'Buoyancy'//time_str(1:MaskSize)
            !     if (gravityProps%type == EQNS_BOD_EXPLICIT) then
            !         call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, txc(1, 1))
            !     else
            !         wrk1d(1:jmax, 1) = 0.0_wp
            !         call Gravity_Source(gravityProps, imax, jmax, kmax, s, txc(1, 1), wrk1d)
            !     end if
            !     dummy = 1.0_wp/froude
            !     txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            !     plot_file = 'LaplacianB'//time_str(1:MaskSize)
            !     call OPR_Partial_Z(OPR_P2, imax, jmax, kmax, g(3), txc(1, 1), txc(1, 4), txc(1, 5))
            !     call OPR_Partial_Y(OPR_P2, imax, jmax, kmax, g(2), txc(1, 1), txc(1, 3), txc(1, 5))
            !     call OPR_Partial_X(OPR_P2, imax, jmax, kmax, g(1), txc(1, 1), txc(1, 2), txc(1, 5))
            !     txc(1:isize_field, 2) = txc(1:isize_field, 2) + txc(1:isize_field, 3) + txc(1:isize_field, 4)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

            !     plot_file = 'GradientRi'//time_str(1:MaskSize)
            !     call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 1), txc(1, 2))
            !     call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), q(1, 1), txc(1, 3))
            !     txc(1:isize_field, 2) = abs(txc(1:isize_field, 2))/(txc(1:isize_field, 3)**2.0 + small_wp)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

            !     plot_file = 'Pressure'//time_str(1:MaskSize)
            !     bbackground = 0.0_wp
            !     call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), DCMP_TOTAL)
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 1))

            !     plot_file = 'PressureGradientY'//time_str(1:MaskSize)
            !     call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), txc(1, 1), txc(1, 2))
            !     call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 2))

            ! end if

            ! ! ###################################################################
            ! ! Total stress tensor
            ! ! ###################################################################
            ! if (opt_vec(iv) == iscal_offset + 20) then ! Total stress tensor
            !     call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 7), txc(1, 1), txc(1, 2), txc(1, 3), DCMP_TOTAL) ! pressure in txc(1,7)
            !     call VISUALS_ACCUMULATE_FIELDS(q, txc(1, 7), txc(1, 8), txc(1, 6))            ! avg vel. + pre. in time
            !     if (it == itime_size) then
            !         call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
            !                               txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            !         txc(1:isize_field, 1) = 2.0_wp*visc*txc(1:isize_field, 1) - txc(1:isize_field, 7)
            !         txc(1:isize_field, 2) = 2.0_wp*visc*txc(1:isize_field, 2) - txc(1:isize_field, 7)
            !         txc(1:isize_field, 3) = 2.0_wp*visc*txc(1:isize_field, 3) - txc(1:isize_field, 7)
            !         txc(1:isize_field, 4) = 2.0_wp*visc*txc(1:isize_field, 4)
            !         txc(1:isize_field, 5) = 2.0_wp*visc*txc(1:isize_field, 5)
            !         txc(1:isize_field, 6) = 2.0_wp*visc*txc(1:isize_field, 6)
            !         if (imode_ibm == 1) then
            !             txc(:, 8) = eps
            !             plot_file = 'EpsSolid'//time_str(1:MaskSize)
            !             call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 1, subdomain, txc(1, 8))
            !         end if
            !         plot_file = 'StressTensor'
            !         if (itime_size > 1) plot_file = trim(adjustl(plot_file))//'Avg'
            !         plot_file = trim(adjustl(plot_file))//time_str(1:MaskSize)
            !         call IO_WRITE_VISUALS(plot_file, imax, jmax, kmax, 6, subdomain, txc(1, 1))
            !     end if

            ! end if

        end do

    end do

    call TLab_Stop(0)

contains
! !########################################################################
! !# Accumulate itime_size fields before processing the temp. avg. fields
! !########################################################################
!     subroutine VISUALS_ACCUMULATE_FIELDS(q, p, q_avg, p_avg)

!         real(wp), dimension(isize_field, 3), intent(inout) :: q     ! inst. vel. fields
!         real(wp), dimension(isize_field), intent(inout) :: p     ! inst. pre. fields
!         real(wp), dimension(isize_field, 3), intent(inout) :: q_avg ! time avg. vel. fields
!         real(wp), dimension(isize_field), intent(inout) :: p_avg ! time avg. pre. fields

!         integer(wi) :: is
!         ! ================================================================== !
!         ! if only one iteration is chosen, do nothing
!         if (itime_size > 1) then
!             if (it == 1) then
!                 q_avg(:, :) = q(:, :); p_avg(:) = p(:)
!             else if (it > 1) then
!                 do is = 1, 3
!                     q_avg(:, is) = q_avg(:, is) + q(:, is)
!                 end do
!                 p_avg(:) = p_avg(:) + p(:)
!             end if
!             if (it == itime_size) then
!                 q = q_avg/itime_size
!                 p = p_avg/itime_size
!             end if
!         end if

!         return
!     end subroutine VISUALS_ACCUMULATE_FIELDS

!########################################################################
!########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

    subroutine IO_WRITE_VISUALS(fname, nx, ny, nz, nfield, subdomain, field)
        character(len=*) fname
        integer(wi) nx, ny, nz, nfield, subdomain(6)
        real(wp), intent(inout) :: field(isize_txc_field, nfield)

        ! -------------------------------------------------------------------
        integer(wi) sizes(5), nx_aux, ny_aux, nz_aux, ifield, i
#ifdef USE_MPI
        character*32 varname(16)
#endif
        character*32 name
        integer(wi) iflag_mode
        integer, parameter :: i1 = 1

        ! ###################################################################
        sizes(5) = nfield

        nx_aux = subdomain(2) - subdomain(1) + 1
        ny_aux = subdomain(4) - subdomain(3) + 1
        nz_aux = subdomain(6) - subdomain(5) + 1

        iflag_mode = 0 ! default
        sizes(1) = isize_txc_field     ! array size
        sizes(2) = 1                   ! lower bound
        ! if (subdomain(2) - subdomain(1) + 1 == g(1)%size .and. &
        !     subdomain(6) - subdomain(5) + 1 == 1) then! xOy plane
        !     iflag_mode = IO_SUBARRAY_VISUALS_XOY
        !     sizes(3) = ny_aux*nx     ! upper bound
        !     sizes(4) = 1              ! stride

        ! else if (subdomain(6) - subdomain(5) + 1 == g(3)%size .and. &
        !          subdomain(2) - subdomain(1) + 1 == 1) then! zOy plane
        !     iflag_mode = IO_SUBARRAY_VISUALS_YOZ
        !     sizes(3) = ny_aux*nx*nz ! upper bound
        !     sizes(4) = nx             ! stride

        ! else
        if (subdomain(2) - subdomain(1) + 1 == x%size .and. &
            subdomain(4) - subdomain(3) + 1 == y%size) then! xOy blocks
            iflag_mode = IO_SUBARRAY_VISUALS_XOY
            sizes(3) = nx*ny*nz_aux ! upper bound
            sizes(4) = 1              ! stride

        end if

        ! ###################################################################
        select case (opt_format)
        case (FORMAT_GENERAL)
            if (nfield > 1 .and. isize_txc_field > nx*ny*nz) then ! IO_Write_Fields expects field to be aligned by nx*ny*nz (instead of isize_txc_field)
                do ifield = 2, nfield
                    do i = 1, nx*ny*nz
                        field((ifield - 1)*nx*ny*nz + i, 1) = field(i, ifield)
                    end do
                end do
            end if
            call IO_Write_Fields(fname, nx, ny, nz, itime, nfield, field)

        case (FORMAT_SINGLE)
#ifdef USE_MPI
            if (nz_aux /= nz) then
                do ifield = 1, nfield
                    call REDUCE_BLOCK_INPLACE(nx, ny, nz, i1, i1, subdomain(5), i1, nx, ny, nz_aux, field(1, ifield), wrk1d)
                end do
            end if

            varname = ''
            if (nfield > 1) then
                do ifield = 1, nfield; write (varname(ifield), *) ifield; varname(ifield) = trim(adjustl(varname(ifield)))
                end do
            end if
            call IO_Write_Subarray(io_subarrays(iflag_mode), fname, varname, field, sizes)
#else
            do ifield = 1, nfield
                call REDUCE_BLOCK_INPLACE(nx, ny, nz, subdomain(1), subdomain(3), subdomain(5), nx_aux, ny_aux, nz_aux, field(1, ifield), wrk1d)

                if (nfield > 1) then
                    write (name, *) ifield; name = trim(adjustl(fname))//'.'//trim(adjustl(name))
                else
                    name = fname
                end if
#include "tlab_open_file.h"
                write (LOC_UNIT_ID) SNGL(field(1:nx_aux*ny_aux*nz_aux, ifield))
                close (LOC_UNIT_ID)

            end do

#endif

        end select

        return
    end subroutine IO_WRITE_VISUALS

    ! ###################################################################
    ! ###################################################################
#ifdef USE_MPI

    subroutine Create_IO_Subarrays()
        ! -----------------------------------------------------------------------
        integer(wi) ny_loc

        ! #######################################################################
        io_subarrays(:)%active = .false.
        io_subarrays(:)%offset = 0
        io_subarrays(:)%precision = IO_TYPE_SINGLE

        ny_loc = subdomain(6) - subdomain(5) + 1

        ! ###################################################################
        ! ! Saving full vertical xOz planes; using subdomain(3) to define the plane
        ! if (ims_pro_j == ((subdomain(3) - 1)/jmax)) io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%active = .true.
        ! io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%communicator = ims_comm_x
        ! io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%subarray = IO_Create_Subarray_XOZ(imax, ny_loc, MPI_REAL4)

        ! ! Saving full vertical yOz planes; using subiddomain(1) to define the plane
        ! if (ims_pro_i == ((subdomain(1) - 1)/imax)) io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%active = .true.
        ! io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%communicator = ims_comm_y
        ! io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%subarray = IO_Create_Subarray_YOZ(ny_loc, kmax, MPI_REAL4)

        ! Saving full blocks xOy planes
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%communicator = MPI_COMM_WORLD
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%subarray = IO_Create_Subarray_XOY(imax, ny_loc, kmax, MPI_REAL4)

        return
    end subroutine Create_IO_Subarrays

#endif

end program VISUALS
