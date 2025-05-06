subroutine OPR_Check()
    use TLab_Constants, only: lfile, wp, wi, fmt_r
    use TLab_Memory, only: isize_field, isize_txc_field, inb_flow_array, inb_txc
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Arrays
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
#ifdef USE_MPI
    use TLab_Time, only: itime
    use mpi_f08
    use TLabMPI_VARS, only: ims_err
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
#endif
    use TLab_Grid, only: x, y, z
    use OPR_Fourier

    implicit none

    ! -------------------------------------------------------------------
    real(wp) residual
    integer(wi) t_srt, t_end, t_dif, PROC_CYCLES, MAX_CYCLES
    real(wp) norm
    character*64 str
    character*256 line

#ifdef USE_MPI
    real(wp) dummy
    integer(wi) idummy
#endif

    complex(wp), pointer :: c_tmp1(:) => null(), c_tmp2(:) => null()

    ! ###################################################################
    if (inb_flow_array < 3 .or. inb_txc < 2) return ! not enough memory

    ! Create random array
    call random_number(q(1:isize_field, 1))

    ! ###################################################################
    ! MPI transposition
#ifdef USE_MPI
    if (ims_npro_i > 1) then            ! MPI transposition along OX
        call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
        call TLabMPI_Trp_ExecI_Forward(q(:, 1), wrk3d, tmpi_plan_dx)
        call TLabMPI_Trp_ExecI_Backward(wrk3d, q(:, 2), tmpi_plan_dx)
        call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        write (str, fmt_r) real(t_dif, wp)/PROC_CYCLES

        dummy = maxval(abs(q(1:isize_field, 1) - q(1:isize_field, 2)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

        write (line, fmt_r) residual
        line = 'Checking MPI-x transposition. Residual ' &
               //trim(adjustl(line))//'. Max. elapsed time '//trim(adjustl(str))//' sec.'
        call TLab_Write_ASCII(lfile, line)

    end if

    if (ims_npro_j > 1) then            ! MPI transposition along Oy
        call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
        idummy = itime; itime = -1      ! set itime to -1 for this call to trigger interruption
        call TLabMPI_Trp_ExecJ_Forward(q(:, 1), wrk3d, tmpi_plan_dy)
        itime = idummy
        call TLabMPI_Trp_ExecJ_Backward(wrk3d, q(:, 2), tmpi_plan_dy)
        call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        write (str, fmt_r) real(t_dif, wp)/PROC_CYCLES

        dummy = maxval(abs(q(1:isize_field, 1) - q(1:isize_field, 2)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

        write (line, fmt_r) residual
        line = 'Checking MPI-y transposition. Residual ' &
               //trim(adjustl(line))//'. Max. elapsed time '//trim(adjustl(str))//' sec.'
        call TLab_Write_ASCII(lfile, line)

    end if
#endif

    ! ###################################################################
    ! FFT
    call c_f_pointer(c_loc(txc(:, 1)), c_tmp1, shape=[isize_txc_field/2])
    call c_f_pointer(c_loc(txc(:, 2)), c_tmp2, shape=[isize_txc_field/2])

    call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
    ! call OPR_Fourier_XY_Forward(q(:, 1), c_tmp1, c_tmp2)
    ! call OPR_Fourier_XY_Backward(c_tmp1, q(:, 2), c_tmp2)
    call OPR_Fourier_Forward(q(:, 1), c_tmp1, c_tmp2)
    call OPR_Fourier_Backward(c_tmp1, q(:, 2), c_tmp2)
    call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)

#ifdef USE_MPI
    idummy = t_end - t_srt
    call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
#else
    t_dif = t_end - t_srt
#endif
    write (str, fmt_r) real(t_dif, wp)/PROC_CYCLES

    if (fft_y_on) then
        ! norm = 1.0_wp/real(x%size*y%size, wp)
        norm = 1.0_wp/real(x%size*y%size*z%size, wp)
    else
        ! norm = 1.0_wp/real(x%size, wp)
        norm = 1.0_wp/real(x%size*z%size, wp)
    end if

#ifdef USE_MPI
    dummy = maxval(abs(norm*q(1:isize_field, 2) - q(1:isize_field, 1)))
    call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
#else
    residual = maxval(abs(norm*q(1:isize_field, 2) - q(1:isize_field, 1)))
#endif

    write (line, fmt_r) residual
    line = 'Checking FFT routines. Residual ' &
           //trim(adjustl(line))//'. Max. elapsed time '//trim(adjustl(str))//' sec.'
    call TLab_Write_ASCII(lfile, line)

    return
end subroutine OPR_Check
