#include "tlab_error.h"

module Statistics
    use TLab_Constants, only: MAX_AVG_TEMPORAL, wp, wi, small_wp
    use TLab_Memory, only: kmax
    ! use Thermodynamics
    implicit none
    private

    public :: Statistics_Initialize, Statistics_Compute

    ! -------------------------------------------------------------------
    logical :: stats_averages = .false.
    logical :: stats_pdfs = .false.
    logical :: stats_intermittency = .false.
    real(wp), allocatable :: mean(:, :)

contains

    ! ###################################################################
    ! ###################################################################
    subroutine Statistics_Initialize(inifile)
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start

        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        ! ###################################################################
        ! Read
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Statistics'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Averages=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#Pdfs=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#Intermittency=<yes/no>')

        call ScanFile_Char(bakfile, inifile, block, 'Averages', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') stats_averages = .true.

        call ScanFile_Char(bakfile, inifile, block, 'Pdfs', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') stats_pdfs = .true.

        call ScanFile_Char(bakfile, inifile, block, 'Intermittency', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') stats_intermittency = .true.

        ! ###################################################################
        ! Memory management
        call TLab_Allocate_Real(__FILE__, mean, [kmax, MAX_AVG_TEMPORAL], 'mean')

        return
    end subroutine Statistics_Initialize

    !########################################################################
    !########################################################################
    subroutine Statistics_Compute()
        ! use TLab_Pointers, only: pointers_dt
        ! use FDM, only: g
        use TLab_Memory, only: inb_scal_array, isize_field
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
        use TLab_WorkFlow, only: scal_on
        use TLab_Arrays
        ! use NavierStokes, only: froude, schmidt
        ! use TLab_Time, only: itime !, rtime
        ! use Thermo_Anelastic
        use NSE_Pressure
        use DNS_ARRAYS
        ! use PARTICLE_VARS
        ! use PARTICLE_ARRAYS
        use FI_VORTICITY_EQN

        ! -------------------------------------------------------------------
        ! real(wp) dummy, amin(16), amax(16)
        integer is !, idummy, nbins, ibc(16), nfield
        ! integer(wi) ij
        ! type(pointers_dt) vars(16)
        ! character*32 fname !, gatename(1)
        ! character*64 str
        ! integer(1) igate
        ! integer(1), allocatable, save :: gate(:)

        ! ###################################################################
        ! Calculate pressure
        if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            ! call NSE_Pressure_Incompressible(q, s, txc(1, 3), txc(1, 1), txc(1, 2), txc(1, 4))
            call NSE_Pressure_Incompressible(q, s, txc(:, 3), txc(:, 4), txc(:, 1), txc(:, 2))
        end if

        ! ###################################################################
        ! Intermittency
        ! ###################################################################
        ! if (stats_intermittency) then
        !     allocate (gate(isize_field))
        !     call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 4))

        !     ! calculate vorticity gate based on 0.1% threshold
        !     call MINMAX(imax, jmax, kmax, txc(1, 1), amin(1), amax(1))
        !     amin(1) = 1.0e-3_wp*1.0e-3_wp*amax(1)
        !     do ij = 1, isize_field
        !         if (txc(ij, 1) > amin(1)) then; gate(ij) = 1  ! gate array
        !         else; gate(ij) = 0
        !         end if
        !     end do
        !     nfield = 1; gatename(1) = 'Vorticity'

        !     write (fname, *) itime; fname = 'int'//trim(adjustl(fname))
        !     call INTER_N_XZ(fname, itime, rtime, imax, jmax, kmax, nfield, gatename, gate, g(2)%nodes, mean)

        !     deallocate (gate)
        ! end if

        ! ###################################################################
        ! Unconditional plane PDFs
        ! ###################################################################
        ! if (stats_pdfs) then
        !     nfield = 0
        !     nfield = nfield + 1; vars(nfield)%field => q(:, 1); vars(nfield)%tag = 'u'
        !     nfield = nfield + 1; vars(nfield)%field => q(:, 2); vars(nfield)%tag = 'v'
        !     nfield = nfield + 1; vars(nfield)%field => q(:, 3); vars(nfield)%tag = 'w'
        !     if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        !         nfield = nfield + 1; vars(nfield)%field => txc(:, 3); vars(nfield)%tag = 'p'
        !     else
        !         nfield = nfield + 1; vars(nfield)%field => q(:, 6); vars(nfield)%tag = 'p'
        !         nfield = nfield + 1; vars(nfield)%field => q(:, 5); vars(nfield)%tag = 'r'
        !         nfield = nfield + 1; vars(nfield)%field => q(:, 7); vars(nfield)%tag = 't'
        !     end if

        !     do is = 1, inb_scal_array
        !         nfield = nfield + 1; vars(nfield)%field => s(:, is); vars(nfield)%tag = 's'
        !         write (str, *) is; vars(nfield)%tag = trim(adjustl(vars(nfield)%tag))//trim(adjustl(str))
        !     end do

        !     ibc(1:nfield) = 2 ! BCs in the calculation of the PDFs
        !     igate = 0         ! no intermittency partition

        !     nbins = 32
        !     write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))
        !     call PDF1V_N(fname, rtime, imax, jmax, kmax, &
        !                  nfield, nbins, ibc, amin, amax, vars, igate, wrk3d, g(2)%nodes, txc)

        ! end if

        ! ###################################################################
        ! Plane averages
        ! ###################################################################
        if (stats_averages) then
            if (scal_on) then
                do is = 1, inb_scal_array          ! All, prognostic and diagnostic fields in array s
                    hq(1:isize_field, 3) = txc(1:isize_field, 3) ! Pass the pressure in hq3
                    call AVG_SCAL_XZ(is, q, s, s(1, is), &
                                     txc(1, 1), txc(1, 2), txc(1, 4), txc(1, 5), txc(1, 6), hq(1, 3), mean)
                end do

            end if

            call AVG_FLOW_XZ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), hq(1, 1), hq(1, 2), hq(1, 3), &
                             mean)
        end if

        return
    end subroutine Statistics_Compute

end module Statistics
