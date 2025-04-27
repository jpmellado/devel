#include "dns_error.h"

module OPR_Fourier
    use TLab_Constants, only: wp, wi, pi_wp, efile
    use TLab_Memory, only: isize_txc_field, isize_txc_dimz
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_Arrays, only: wrk3d
    use TLab_Pointers_C, only: c_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Grid
    use, intrinsic :: iso_c_binding
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
    use TLabMPI_Transpose
#endif
    implicit none
    private

    public :: OPR_Fourier_Initialize
    public :: OPR_Fourier_X_Forward, OPR_Fourier_X_Backward     ! Main routines, used in Poisson solver; optimize
    public :: OPR_Fourier_Y_Forward, OPR_Fourier_Y_Backward
    public :: OPR_Fourier_SetPSD

    logical, public :: fft_x_on = .false.                       ! Switch on and off the FFT in the corresponding directions.
    logical, public :: fft_y_on = .false.
    logical, public :: fft_z_on = .false.

    ! -----------------------------------------------------------------------
    integer(wi) size_fft_x, size_fft_y, size_fft_z

    type(c_ptr) :: fft_plan_fx, fft_plan_bx
    type(c_ptr) :: fft_plan_fy, fft_plan_by!, fft_plan_fy1d, fft_plan_by1d
    type(c_ptr) :: fft_plan_fz, fft_plan_bz
#ifdef USE_MPI
    type(tmpi_transpose_dt), public :: tmpi_plan_fftx
    type(tmpi_transpose_dt), public :: tmpi_plan_ffty
#endif

    complex(wp), pointer :: c_out(:) => null()
    real(wp), pointer :: r_out(:) => null()

contains
    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_Initialize()
#ifdef USE_MPI
        use mpi_f08, only: MPI_DOUBLE_COMPLEX
#endif
        use TLab_Arrays, only: txc

#ifdef USE_FFTW
#include "fftw3.f03"
#endif
        ! -----------------------------------------------------------------------
        integer(wi) stride, offset, nlines

        ! #######################################################################
#ifndef USE_FFTW
        call TLab_Write_ASCII(efile, __FILE__//'. FFTW needed for POISSON solver.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
#endif

        if (mod(imax, 2) /= 0) then
            call TLab_Write_ASCII(efile, __FILE__//'. Imax must be a multiple of 2 for the FFT operations.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! -----------------------------------------------------------------------
        ! Ox direction
        size_fft_x = size(x(:))
        if (size_fft_x > 1) then
            fft_x_on = .true.

            nlines = jmax*kmax
            offset = size_fft_x/2 + 1

#ifdef USE_MPI
            if (ims_npro_i > 1) then
                ! Extended with the Nyquist frequency
                tmpi_plan_fftx = TLabMPI_Trp_PlanI(imax/2 + 1, nlines, &
                                                   locType=MPI_DOUBLE_COMPLEX, &
                                                   message='extended Ox FFTW in Poisson solver.')

                nlines = tmpi_plan_dx%nlines
                offset = (imax/2 + 1)*ims_npro_i

            end if
#endif

#ifdef _DEBUG
#else
            call dfftw_plan_many_dft_r2c(fft_plan_fx, 1, size_fft_x, nlines, &
                                         txc(:, 1), size_fft_x, 1, size_fft_x, &
                                         wrk3d, size_fft_x/2 + 1, 1, offset, &
                                         FFTW_MEASURE)
            call dfftw_plan_many_dft_c2r(fft_plan_bx, 1, size_fft_x, nlines, &
                                         txc(:, 1), size_fft_x/2 + 1, 1, offset, &
                                         wrk3d, size_fft_x, 1, size_fft_x, &
                                         FFTW_MEASURE)
#endif

        end if

        ! -----------------------------------------------------------------------
        ! Oy direction
        size_fft_y = size(y(:))

        if (size_fft_y > 1) then
            fft_y_on = .true.

            nlines = (imax/2 + 1)*kmax

#ifdef USE_MPI
            if (ims_npro_j > 1) then
                tmpi_plan_ffty = TLabMPI_Trp_PlanJ(jmax, nlines, &
                                                   locType=MPI_DOUBLE_COMPLEX, &
                                                   message='Oz FFTW in Poisson solver.')
                nlines = tmpi_plan_fftz%nlines

            end if
#endif

            stride = nlines

#ifdef _DEBUG
#else
            call dfftw_plan_many_dft(fft_plan_fy, 1, size_fft_y, nlines, &
                                     txc(:, 1), size_fft_y, stride, 1, &
                                     wrk3d, size_fft_y, stride, 1, &
                                     FFTW_FORWARD, FFTW_MEASURE)

            call dfftw_plan_many_dft(fft_plan_by, 1, size_fft_y, nlines, &
                                     txc(:, 1), size_fft_y, stride, 1, &
                                     wrk3d, size_fft_y, stride, 1, &
                                     FFTW_BACKWARD, FFTW_MEASURE)
#endif
        end if

        ! -----------------------------------------------------------------------
        ! Oz direction
        size_fft_z = size(z(:))
        if (size_fft_z > 1) then
            fft_z_on = .true.

            nlines = (imax/2 + 1)*jmax
            stride = (imax/2 + 1)*jmax

#ifdef _DEBUG
#else
            call dfftw_plan_many_dft(fft_plan_fz, 1, size_fft_z, nlines, &
                                     txc(:, 1), size_fft_z, stride, 1, &
                                     wrk3d, size_fft_z, stride, 1, &
                                     FFTW_FORWARD, FFTW_MEASURE)

            call dfftw_plan_many_dft(fft_plan_bz, 1, size_fft_z, nlines, &
                                     txc(:, 1), size_fft_z, stride, 1, &
                                     wrk3d, size_fft_z, stride, 1, &
                                     FFTW_BACKWARD, FFTW_MEASURE)
#endif
        end if

        return
    end subroutine OPR_Fourier_Initialize

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_X_Forward(in, out)
        real(wp), intent(in) :: in(*)
        complex(wp), intent(out) :: out(*)

        target out

        ! -----------------------------------------------------------------------

        ! #######################################################################

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call c_f_pointer(c_loc(out), r_out, shape=[isize_txc_field])

            call TLabMPI_Trp_ExecI_Forward(in(:), r_out(:), tmpi_plan_dx)

            call dfftw_execute_dft_r2c(fft_plan_fx, r_out(:), c_wrk3d)

            call TLabMPI_Trp_ExecI_Backward(c_wrk3d, out(:), tmpi_plan_fftx)

            nullify (r_out)

        else
#endif
            call dfftw_execute_dft_r2c(fft_plan_fx, in, out)

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Fourier_X_Forward

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_X_Backward(in, out)
        complex(wp), intent(in) :: in(*)
        real(wp), intent(out) :: out(*)

        ! -----------------------------------------------------------------------
        !########################################################################
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call c_f_pointer(c_loc(out), c_out, shape=[isize_txc_field])

            call TLabMPI_Trp_ExecI_Forward(in(:), c_out(:, 1), tmpi_plan_fftx)

            call dfftw_execute_dft_c2r(fft_plan_bx, c_out, wrk3d)

            call TLabMPI_Trp_ExecI_Backward(wrk3d, out(:), tmpi_plan_dx) !tmpi_plan_fftx1)

            nullify (r_in, c_out)

        else
#endif
            call dfftw_execute_dft_c2r(fft_plan_bx, in, out)

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Fourier_X_Backward

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_Z_Forward(in, out)
        complex(wp), intent(inout), target :: in(*)     ! in might be overwritten
        complex(wp), intent(out), target :: out(*)

        ! -----------------------------------------------------------------------

        ! #######################################################################
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_Trp_ExecK_Forward(in, out, tmpi_plan_fftz)

            call dfftw_execute_dft(fft_plan_fz, out, c_wrk3d)

            call TLabMPI_Trp_ExecK_Backward(c_wrk3d, out, tmpi_plan_fftz)

        else
#endif
            call dfftw_execute_dft(fft_plan_fz, in, out)

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Fourier_Z_Forward

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_Z_Backward(in, out)
        complex(wp), intent(inout) :: in(*)
        complex(wp), intent(out) :: out(*)

        ! #######################################################################
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_Trp_ExecK_Forward(in, out, tmpi_plan_fftz)

            call dfftw_execute_dft(fft_plan_bz, out, c_wrk3d)

            call TLabMPI_Trp_ExecK_Backward(c_wrk3d, out, tmpi_plan_fftz)

        else
#endif
            call dfftw_execute_dft(fft_plan_bz, in, out)
#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Fourier_Z_Backward

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_SetPSD(nx, ny, nz, u, psd, locPhase)
        use Distributions
        use FDM, only: g

        integer(wi) nx, ny, nz
        complex(wp), intent(inout) :: u(nx/2 + 1, ny, nz)
        type(distributions_dt), intent(in) :: psd
        real(wp), intent(in), optional :: locPhase(nx/2 + 1, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) i, j, iglobal, jglobal, kglobal
        real(wp) pow_dst, pow_org, phase
        real(wp) f, fi, fj, fk

        ! #######################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= size_fft_z/2 + 1) then
                fk = real(kglobal - 1, wp)/g(3)%scale
            else
                fk = -real(size_fft_z + 1 - kglobal, wp)/g(3)%scale
            end if

            do j = 1, ny
                jglobal = j ! No MPI decomposition along Oy
                if (jglobal <= size_fft_y/2 + 1) then
                    fj = real(jglobal - 1, wp)/g(2)%scale
                else
                    fj = -real(size_fft_y + 1 - jglobal, wp)/g(2)%scale
                end if

                do i = 1, nx/2 + 1
#ifdef USE_MPI
                    iglobal = i + ims_offset_i/2
#else
                    iglobal = i
#endif
                    fi = real(iglobal - 1, wp)/g(1)%scale

                    f = sqrt(fi**2 + fj**2 + fk**2)

                    ! -----------------------------------------------------------------------
                    ! target psd
                    pow_dst = Distributions_Compute(psd, f)

                    if (f == 0.0_wp) then
                        pow_dst = 0.0_wp

                    else
                        if (size_fft_y == 1 .or. size_fft_z == 1) then ! 2D spectrum
                            pow_dst = pow_dst/(pi_wp*f)
                        else
                            pow_dst = pow_dst/(2*pi_wp*f**2)
                        end if

                    end if

                    ! -----------------------------------------------------------------------
                    ! phase and scaling of complex data
                    if (pow_dst > 0.0_wp) pow_dst = sqrt(pow_dst)

                    if (present(locPhase)) then
                        if (iglobal == 1 .or. iglobal == size_fft_x/2 + 1) then
                            phase = 0.0_wp
                        else
                            phase = (locPhase(i, j, k) - 0.5_wp)*2.0_wp*pi_wp
                        end if
                        u(i, j, k) = cmplx(pow_dst*cos(phase), pow_dst*sin(phase), wp)

                    else
                        pow_org = abs(u(i, j, k))

                        if (pow_org > 0.0_wp) pow_dst = pow_dst/pow_org

                        u(i, j, k) = u(i, j, k)*pow_dst

                    end if

                end do
            end do
        end do

        return
    end subroutine OPR_Fourier_SetPSD

end module OPR_Fourier
