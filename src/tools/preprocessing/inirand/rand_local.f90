#include "tlab_error.h"

module RAND_LOCAL
    use TLab_Constants, only: wp, wi, big_wp
    use TLab_Constants, only: efile, lfile
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_txc_field, inb_txc
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif
    use TLab_Grid, only: x, y, z
    use Averages, only: AVG1V2D
    use Distributions
    use OPR_Fourier
    implicit none

    type(distributions_dt) :: psd
    type(distributions_dt) :: pdf
    real(wp) :: ucov(6)
    integer(wi) :: seed          ! Random number generator

    ! -------------------------------------------------------------------
    integer(wi) i

contains

    ! ###################################################################
    subroutine Inirand_Initialize_Parameters(inifile)
        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer(wi) :: idummy
        real(wp) :: rdummy(6)

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        
        block = 'Broadband'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Spectrum=<none/uniform/quartic/quadratic/gaussian>')
        call TLab_Write_ASCII(bakfile, '#f0=<frequencies>')
        call TLab_Write_ASCII(bakfile, '#Distribution=<none/uniform/gaussian>')
        call TLab_Write_ASCII(bakfile, '#Covariance=<Rxx,Ryy,Rzz,Rxy,Rxz,Ryz>')
        call TLab_Write_ASCII(bakfile, '#Seed=<random seed>')

        call ScanFile_Int(bakfile, inifile, block, 'Seed', '7', seed)
#ifdef USE_MPI
        seed = seed + ims_pro         ! seed for random generator
#endif
        seed = -abs(seed)

        call ScanFile_Char(bakfile, inifile, block, 'Spectrum', 'quartic', sRes)
        if (trim(adjustl(sRes)) == 'none') then; psd%type = TYPE_DF_NONE
        else if (trim(adjustl(sRes)) == 'uniform') then; psd%type = TYPE_DF_UNIFORM
        else if (trim(adjustl(sRes)) == 'quartic') then; psd%type = TYPE_DF_QUARTIC
        else if (trim(adjustl(sRes)) == 'quadratic') then; psd%type = TYPE_DF_QUADRATIC
        else if (trim(adjustl(sRes)) == 'gaussian') then; psd%type = TYPE_DF_GAUSSIAN
        end if

        psd%parameters(1) = 0.0_wp
        psd%parameters(2:3) = [0.0_wp, big_wp]
        call ScanFile_Char(bakfile, inifile, block, 'f0', '1.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, psd%parameters)
        psd%mean = psd%parameters(1)
        psd%parameters(1) = psd%parameters(2)
        psd%parameters(2) = psd%parameters(3)

        call ScanFile_Real(bakfile, inifile, block, 'Sigma', '-1.0', psd%sigma)
        if (psd%sigma < 0.0_wp) psd%sigma = psd%mean/6.0_wp    ! Default value

        call ScanFile_Char(bakfile, inifile, block, 'Distribution', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; pdf%type = TYPE_DF_NONE
        else if (trim(adjustl(sRes)) == 'uniform') then; pdf%type = TYPE_DF_UNIFORM
        else if (trim(adjustl(sRes)) == 'gaussian') then; pdf%type = TYPE_DF_GAUSSIAN
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Distribution type unknown.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ucov(1:3) = 1.0_wp ! diagonal terms
        ucov(4:6) = 0.0_wp ! off-diagonal terms
        call ScanFile_Char(bakfile, inifile, block, 'Covariance', '-1', sRes)
        if (trim(adjustl(sRes)) /= '-1') then
            idummy = 6
            call LIST_REAL(sRes, idummy, rdummy)
            if (idummy == 6) then; ucov(1:6) = rdummy(1:6)
            else
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Incorrect number of variances.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if
! CALL TLab_Write_ASCII(bakfile,'Covariance matrix:')
! WRITE(sRes,'(3E11.4)') ucov(1), ucov(4), ucov(5); CALL TLab_Write_ASCII(bakfile,sRes)
! WRITE(sRes,'(3E11.4)') ucov(4), ucov(2), ucov(6); CALL TLab_Write_ASCII(bakfile,sRes)
! WRITE(sRes,'(3E11.4)') ucov(5), ucov(6), ucov(3); CALL TLab_Write_ASCII(bakfile,sRes)

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        inb_txc = 3

        return
    end subroutine Inirand_Initialize_Parameters

    ! ###################################################################
    subroutine RAND_FIELD(variance, a, tmp1, tmp2, tmp3)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

        real(wp), intent(IN) :: variance
        real(wp), intent(OUT) :: a(:)
        real(wp), intent(INOUT) :: tmp1(:), tmp2(:), tmp3(:)

        ! -------------------------------------------------------------------
        real(wp) RAN0, RANG
        external RAN0, RANG
        complex(wp), pointer :: c_tmp1(:) => null(), c_tmp3(:) => null()
        target tmp1, tmp3

        ! ###################################################################
        select case (pdf%type)
        case (TYPE_DF_UNIFORM)
            do i = 1, isize_field
                tmp2(i) = RAN0(seed) - 0.5_wp
            end do

        case (TYPE_DF_GAUSSIAN)
            do i = 1, isize_field
                tmp2(i) = RANG(0.0_wp, 1.0_wp, seed)
            end do

        end select

        if (psd%type > 0) then
            call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field/2])
            call c_f_pointer(c_loc(tmp3), c_tmp3, shape=[isize_txc_field/2])
            if (pdf%type > 0) then
                call OPR_Fourier_Forward(tmp2, c_tmp1, c_tmp3)
                call OPR_Fourier_SetPSD(imax, jmax, kmax, c_tmp1, psd)
            else
                do i = 1, isize_txc_field
                    tmp3(i) = RAN0(seed)
                end do
                call OPR_Fourier_SetPSD(imax, jmax, kmax, c_tmp1, psd, locPhase=tmp3)
            end if
            call OPR_Fourier_Backward(c_tmp1, tmp2, c_tmp3)
            nullify (c_tmp1, c_tmp3)
        end if

        call RAND_NORMALIZE(variance, tmp2)
        a(1:isize_field) = tmp2(1:isize_field)

        return
    end subroutine RAND_FIELD

!########################################################################
    subroutine RAND_COVARIANCE(cov, u, v, w)
        real(wp) cov(6)
        real(wp), intent(out) :: u(:), v(:), w(:)

        ! -------------------------------------------------------------------
        real(wp) trace, lambda1, lambda2, alpha, calpha, salpha
        real(wp) rdummy

#define Rxx cov(1)
#define Ryy cov(2)
#define Rzz cov(3)
#define Rxy cov(4)
#define Rxz cov(5)
#define Ryz cov(6)

        ! ###################################################################
        if (y%size > 1) then
            if (Rxy /= 0.0_wp .or. Ryz /= 0.0_wp) then ! only 2D case developed
                call TLab_Write_ASCII(efile, 'Terms Rxy and Ryz not developed yet.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            call RAND_NORMALIZE(Ryy, w)

        end if

        if (Rxz == 0.0_wp) then  ! Diagonal case
            call RAND_NORMALIZE(Rxx, u)
            call RAND_NORMALIZE(Rzz, v)

        else                        ! Nondiagonal case
            ! get eigenvalues
            trace = Rxx + Rzz
            lambda1 = 0.5_wp*(trace + sqrt(trace*trace - 4.0_wp*(Rxx*Rzz - Rxz*Rxz)))
            lambda2 = trace - lambda1

            ! define fields in rotated uncorrelated frame
            call RAND_NORMALIZE(lambda1, u)
            call RAND_NORMALIZE(lambda2, w)

            ! rotate to XZ correlated frame
            alpha = atan((lambda1 - Rxx)/Rxz)
            calpha = cos(alpha)
            salpha = sin(alpha)

            do i = 1, isize_field
                rdummy = calpha*u(i) - salpha*w(i)
                w(i) = salpha*u(i) + calpha*w(i)
                u(i) = rdummy
            end do

        end if

        return

#undef Rxx
#undef Ryy
#undef Rzz
#undef Rxy
#undef Rxz
#undef Ryz

    end subroutine RAND_COVARIANCE

! ###################################################################
    subroutine RAND_NORMALIZE(variance, a)
        real(wp), intent(IN) :: variance
        real(wp), intent(INOUT) :: a(:)

        ! -------------------------------------------------------------------
        real(wp) dummy

        ! ###################################################################
        dummy = AVG1V2D(size(a), 1, 1, 1, 1, a) ! 3D average
        a = a - dummy

        dummy = AVG1V2D(size(a), 1, 1, 1, 2, a) ! 3D average
        if (dummy > 0.0_wp) then
            dummy = sqrt(variance/dummy)
            a = a*dummy
        end if

        return
    end subroutine RAND_NORMALIZE

end module RAND_LOCAL
