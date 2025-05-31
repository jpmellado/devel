#include "tlab_error.h"

module LargeScaleForcing
    use TLab_Constants, only: wp, wi, efile, MAX_PROF, MAX_VARS, MAX_PARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Grid, only: z
    implicit none
    private

    public :: LargeScaleForcing_Initialize
    public :: LargeScaleForcing_Subsidence

    real(wp), allocatable, public :: wbackground(:)
    ! real(wp), allocatable, public :: ubackground(:)
    ! real(wp), allocatable, public :: vbackground(:)

    ! -------------------------------------------------------------------
    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt
    type(term_dt), public, protected :: subsidenceProps

    integer, parameter :: TYPE_SUB_NONE = 0
    integer, parameter, public :: TYPE_SUB_CONSTANT = 1

contains
    !########################################################################
    !########################################################################
    subroutine LargeScaleForcing_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Subsidence'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/ConstantDivergence>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; subsidenceProps%type = TYPE_SUB_NONE
        else if (trim(adjustl(sRes)) == 'constantdivergence') then; subsidenceProps%type = TYPE_SUB_CONSTANT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        subsidenceProps%active = .false.
        if (subsidenceProps%type /= TYPE_SUB_NONE) then
            subsidenceProps%active = .true.

            subsidenceProps%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '0.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, subsidenceProps%parameters)
        end if

        ! -------------------------------------------------------------------
        ! Initialize
        allocate (wbackground(z%size))
        wbackground(:) = z%nodes(:)*subsidenceProps%parameters(1)      ! subsidence velocity

        return
    end subroutine LargeScaleForcing_Initialize

    !########################################################################
    !########################################################################
    subroutine LargeScaleForcing_Subsidence(locProps, nx, ny, nz, a, source)
        use OPR_Partial, only: OPR_Partial_Z, OPR_P1
        use Averages, only: AVG1V2D_V
        use FDM, only: g

        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: a(nx, ny, nz)
        real(wp), intent(out) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) k

        !########################################################################
        select case (locProps%type)
        case (TYPE_SUB_CONSTANT)
            call OPR_Partial_Z(OPR_P1, nx, ny, nz, g(3), a, source)

            do k = 1, nz
                source(:, :, k) = source(:, :, k)*wbackground(k)
            end do

        end select

        return
    end subroutine LargeScaleForcing_Subsidence

end module LargeScaleForcing
