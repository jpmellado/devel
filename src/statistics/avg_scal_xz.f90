#include "tlab_error.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array tmp3 until pressure block
!#
!########################################################################

subroutine AVG_SCAL_XZ(is, q, s, s_local, dsdx, dsdy, dsdz, tmp1, tmp2, tmp3, mean2d)
    use TLab_Constants, only: MAX_AVG_TEMPORAL
    use TLab_Constants, only: efile, lfile, wp, wi
    use TLab_Memory, only: imax, jmax, kmax, inb_flow_array, inb_scal_array
    use TLab_Time, only: itime, rtime
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: p_wrk3d, u, v, w
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Grid, only: z
    use FDM, only: g
    use FDM, only: fdm_Int0
    use OPR_Partial
    use NavierStokes
    use Thermo_Anelastic, only: ribackground, Thermo_Anelastic_Weight_InPlace
    use Thermo_Compressible, only: rho, vis
    use Gravity, only: gravityProps, Gravity_Source
    ! use Rotation, only: coriolis
    use Microphysics
    use Radiation
    use FI_GRADIENT_EQN
    use Averages, only: AVG_IK_V

    implicit none

    integer, intent(in) :: is
    real(wp), intent(in) :: q(imax, jmax, kmax, inb_flow_array)
    real(wp), intent(in) :: s(imax, jmax, kmax, inb_scal_array)
    real(wp), intent(in) :: s_local(imax, jmax, kmax)
    real(wp), intent(inout) :: dsdx(imax, jmax, kmax)
    real(wp), intent(inout) :: dsdy(imax, jmax, kmax)
    real(wp), intent(inout) :: dsdz(imax, jmax, kmax)
    real(wp), intent(inout) :: tmp1(imax, jmax, kmax)
    real(wp), intent(inout) :: tmp2(imax, jmax, kmax)
    real(wp), intent(inout) :: tmp3(imax, jmax, kmax)
    real(wp), intent(inout) :: mean2d(kmax, MAX_AVG_TEMPORAL)

    target q, tmp3

    ! -----------------------------------------------------------------------
    integer, parameter :: MAX_VARS_GROUPS = 10
    integer k, is_loc
    real(wp) diff, c23

    integer ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nv, iv, offset_sources

    character*32 name, groupname(MAX_VARS_GROUPS)
    character*250 line1, varname(MAX_VARS_GROUPS)

    ! Pointers to existing allocated space
    real(wp), dimension(:, :, :), pointer :: p_loc

    integer :: itransport = 0, EQNS_TRANS_SUTHERLAND = 1, EQNS_TRANS_POWERLAW = 2! temporal fix

    ! ###################################################################
    ! Define pointers
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        p_loc => q(:, :, :, 6)
    else
        p_loc => tmp3
    end if

    if (nse_diffusion == EQNS_NONE) then
        diff = 0.0_wp
    else
        diff = visc/schmidt(is)
    end if

    c23 = 2.0_wp/3.0_wp

    ! -----------------------------------------------------------------------
    ! Dependent variables
    ng = 1; ig(ng) = 1
#define rS(k)     mean2d(k,ig(1)  )
#define fS(k)     mean2d(k,ig(1)+1)
#define rS_z(k)   mean2d(k,ig(1)+2)
#define fS_z(k)   mean2d(k,ig(1)+3)
#define rQ(k)     mean2d(k,ig(1)+4)
#define fQ(k)     mean2d(k,ig(1)+5)
    sg(ng) = 6

    groupname(ng) = 'Mean'
    varname(ng) = 'rS fS rS_z fS_z rQ fQ'

    offset_sources = sg(ng)
    if (infraredProps%active(is)) then
        varname(ng) = trim(adjustl(varname(ng)))//' rQrad rFrad rFradUp'
        sg(ng) = sg(ng) + 3
    end if
    if (sedimentationProps%active(is)) then
        varname(ng) = trim(adjustl(varname(ng)))//' rQsed rFsed'
        sg(ng) = sg(ng) + 2
    end if
    if (evaporationProps%active(is)) then
        varname(ng) = trim(adjustl(varname(ng)))//' rQeva'
        sg(ng) = sg(ng) + 1
    end if

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsu(k)    mean2d(k,ig(2)  )
#define Rsv(k)    mean2d(k,ig(2)+1)
#define Rsw(k)    mean2d(k,ig(2)+2)
#define fS2(k)    mean2d(k,ig(2)+3)
#define fS3(k)    mean2d(k,ig(2)+4)
#define fS4(k)    mean2d(k,ig(2)+5)
#define rS2(k)    mean2d(k,ig(2)+6)
#define rS3(k)    mean2d(k,ig(2)+7)
#define rS4(k)    mean2d(k,ig(2)+8)
    sg(ng) = 9

    groupname(ng) = 'Fluctuations'
    varname(ng) = 'Rsu Rsv Rsw fS2 fS3 fS4 rS2 rS3 rS4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rss_t(k)  mean2d(k,ig(3)  )
#define Css(k)    mean2d(k,ig(3)+1)
#define Pss(k)    mean2d(k,ig(3)+2)
#define Ess(k)    mean2d(k,ig(3)+3)
#define Tssz1(k)  mean2d(k,ig(3)+4)
#define Tssz2(k)  mean2d(k,ig(3)+5)
#define Tssz_z(k) mean2d(k,ig(3)+6)
#define Dss(k)    mean2d(k,ig(3)+7)
#define Qss(k)    mean2d(k,ig(3)+8)
    sg(ng) = 9

    groupname(ng) = 'RssBudget'
    varname(ng) = 'Rss_t Css Pss Ess Tssz1 Tssz2 Tssz_z Dss Qss'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsu_t(k)  mean2d(k,ig(4)  )
#define Csu(k)    mean2d(k,ig(4)+1)
#define Psu(k)    mean2d(k,ig(4)+2)
#define Esu(k)    mean2d(k,ig(4)+3)
#define PIsu(k)   mean2d(k,ig(4)+4)
#define Tsuz1(k)  mean2d(k,ig(4)+5)
#define Tsuz2(k)  mean2d(k,ig(4)+6)
#define Tsuz_z(k) mean2d(k,ig(4)+7)
#define Dsu(k)    mean2d(k,ig(4)+8)
#define Gsu(k)    mean2d(k,ig(4)+9)
#define Bsu(k)    mean2d(k,ig(4)+10)
#define Fsu(k)    mean2d(k,ig(4)+11)
#define Qsu(k)    mean2d(k,ig(4)+12)
    sg(ng) = 13

    groupname(ng) = 'RsuBudget'
    varname(ng) = 'Rsu_t Csu Psu Esu PIsu Tsuz1 Tsuz2 Tsuz_z Dsu Gsu Bsu Fsu Qsu'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsv_t(k)  mean2d(k,ig(5)  )
#define Csv(k)    mean2d(k,ig(5)+1)
#define Psv(k)    mean2d(k,ig(5)+2)
#define Esv(k)    mean2d(k,ig(5)+3)
#define PIsv(k)   mean2d(k,ig(5)+4)
#define Tsvz1(k)  mean2d(k,ig(5)+5)
#define Tsvz2(k)  mean2d(k,ig(5)+6)
#define Tsvz_z(k) mean2d(k,ig(5)+7)
#define Dsv(k)    mean2d(k,ig(5)+8)
#define Gsv(k)    mean2d(k,ig(5)+9)
#define Bsv(k)    mean2d(k,ig(5)+10)
#define Fsv(k)    mean2d(k,ig(5)+11)
#define Qsv(k)    mean2d(k,ig(5)+12)
    sg(ng) = 14

    groupname(ng) = 'RsvBudget'
    varname(ng) = 'Rsv_t Csv Psv Esv PIsv Tsvz1 Tsvz2 Tsvz_z Dsv Gsv Bsv Fsv Qsv'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsw_t(k)  mean2d(k,ig(6)  )
#define Csw(k)    mean2d(k,ig(6)+1)
#define Psw(k)    mean2d(k,ig(6)+2)
#define Esw(k)    mean2d(k,ig(6)+3)
#define PIsw(k)   mean2d(k,ig(6)+4)
#define Tswz1(k)  mean2d(k,ig(6)+5)
#define Tswz2(k)  mean2d(k,ig(6)+6)
#define Tswz3(k)  mean2d(k,ig(5)+7)
#define Tswz_z(k) mean2d(k,ig(6)+8)
#define Dsw(k)    mean2d(k,ig(6)+9)
#define Gsw(k)    mean2d(k,ig(6)+10)
#define Bsw(k)    mean2d(k,ig(6)+11)
#define Fsw(k)    mean2d(k,ig(6)+12)
#define Qsw(k)    mean2d(k,ig(6)+13)
    sg(ng) = 13

    groupname(ng) = 'RswBudget'
    varname(ng) = 'Rsw_t Csw Psw Esw PIsw Tswz1 Tswz2 Tsvz3 Tswz_z Dsw Gsw Bsw Fsw Qsw'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define S_x2(k)  mean2d(k,ig(7)  )
#define S_y2(k)  mean2d(k,ig(7)+1)
#define S_z2(k)  mean2d(k,ig(7)+2)
#define S_x3(k)  mean2d(k,ig(7)+3)
#define S_y3(k)  mean2d(k,ig(7)+4)
#define S_z3(k)  mean2d(k,ig(7)+5)
#define S_x4(k)  mean2d(k,ig(7)+6)
#define S_y4(k)  mean2d(k,ig(7)+7)
#define S_z4(k)  mean2d(k,ig(7)+8)
    sg(ng) = 9

    groupname(ng) = 'DerivativeFluctuations'
    varname(ng) = 'S_x2 S_y2 S_z2 S_x3 S_y3 S_z3 S_x4 S_y4 S_z4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
    sg(ng) = 2*inb_scal_array

    groupname(ng) = 'CrossScalars'
    varname(ng) = ''
    do is_loc = 1, inb_scal_array
        write (name, *) is_loc
        varname(ng) = trim(adjustl(varname(ng)))//' Cs'//trim(adjustl(name))
        varname(ng) = trim(adjustl(varname(ng)))//' Css'//trim(adjustl(name))
    end do

    ! -----------------------------------------------------------------------
    ! Auxiliary variables depending on z and t; this last group is not written
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rR(k)       mean2d(k,ig(9))
#define rU(k)       mean2d(k,ig(9)+1)
#define rV(k)       mean2d(k,ig(9)+2)
#define rW(k)       mean2d(k,ig(9)+3)
#define fU(k)       mean2d(k,ig(9)+4)
#define fV(k)       mean2d(k,ig(9)+5)
#define fW(k)       mean2d(k,ig(9)+6)
#define fU_z(k)     mean2d(k,ig(9)+7)
#define fV_z(k)     mean2d(k,ig(9)+8)
#define fW_z(k)     mean2d(k,ig(9)+9)
#define rU_z(k)     mean2d(k,ig(9)+10)
#define rV_z(k)     mean2d(k,ig(9)+11)
#define rW_z(k)     mean2d(k,ig(9)+12)
#define Rwu(k)      mean2d(k,ig(9)+13)
#define Rwv(k)      mean2d(k,ig(9)+14)
#define Rww(k)      mean2d(k,ig(9)+15)
#define Rss_z(k)    mean2d(k,ig(9)+16)
#define Rsu_z(k)    mean2d(k,ig(9)+17)
#define Rsv_z(k)    mean2d(k,ig(9)+18)
#define Rsw_z(k)    mean2d(k,ig(9)+19)
#define Fz(k)       mean2d(k,ig(9)+20)
#define Fz_z(k)     mean2d(k,ig(9)+21)
#define Tau_zz(k)   mean2d(k,ig(9)+22)
#define Tau_zz_z(k) mean2d(k,ig(9)+23)
#define Tau_zx(k)   mean2d(k,ig(9)+24)
#define Tau_zx_z(k) mean2d(k,ig(9)+25)
#define Tau_yz(k)   mean2d(k,ig(9)+26)
#define Tau_yz_z(k) mean2d(k,ig(9)+27)
#define rP(k)       mean2d(k,ig(9)+28)
#define aux(k)      mean2d(k,ig(9)+29)
    sg(ng) = 30

    ! -----------------------------------------------------------------------
    nv = ig(ng) + sg(ng) - 1
    if (MAX_AVG_TEMPORAL < nv) then
        call TLab_Write_ASCII(efile, __FILE__//'. Not enough space in local arrays.')
        call TLab_Stop(DNS_ERROR_AVGTMP)
    end if
    mean2d(:, 1:nv) = 0.0_wp

    ng = ng - 1
    nv = ig(ng) + sg(ng) - 1 ! the last group is not written out

    ! #######################################################################
    write (line1, *) itime; line1 = 'Calculating scal statistics at It'//trim(adjustl(line1))//'...'
    call TLab_Write_ASCII(lfile, line1)

    ! #######################################################################
    ! Preliminary data of velocity and density
    ! #######################################################################
    call AVG_IK_V(imax, jmax, kmax, u, rU(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, v, rV(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, w, rW(1), wrk1d)

    if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        rR(:) = 1.0_wp

        fU(:) = rU(:)
        fV(:) = rV(:)
        fW(:) = rW(:)

    else
        call AVG_IK_V(imax, jmax, kmax, rho, rR(1), wrk1d)

        dsdx = rho*u
        dsdy = rho*v
        dsdz = rho*w
        call AVG_IK_V(imax, jmax, kmax, dsdx, fU(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dsdy, fV(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dsdz, fW(1), wrk1d)
        fU(:) = fU(:)/rR(:)
        fV(:) = fV(:)/rR(:)
        fW(:) = fW(:)/rR(:)

    end if

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rU(1), rU_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rV(1), rV_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rW(1), rW_z(1))

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fU(1), fU_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fV(1), fV_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fW(1), fW_z(1))

    dsdx = w*u
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) dsdx = dsdx*rho
    dsdy = w*v
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) dsdy = dsdy*rho
    dsdz = w*w
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) dsdz = dsdz*rho
    call AVG_IK_V(imax, jmax, kmax, dsdx, Rwu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Rwv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Rww(1), wrk1d)
    Rwu(:) = Rwu(:)/rR(:) - fW(:)*fU(:)
    Rwv(:) = Rwv(:)/rR(:) - fW(:)*fV(:)
    Rww(:) = Rww(:)/rR(:) - fW(:)*fW(:)

    ! #######################################################################
    ! Scalar
    ! #######################################################################
    call AVG_IK_V(imax, jmax, kmax, s_local, rS(1), wrk1d)

    if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        fS(:) = rS(:)
    else
        p_wrk3d = rho*s_local
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, fS(1), wrk1d)
        fS(:) = fS(:)/rR(:)
    end if

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rS(1), rS_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fS(1), fS_z(1))

    ! -----------------------------------------------------------------------
    ! Moments
    do k = 1, kmax
        p_wrk3d(:, :, k) = s_local(:, :, k) - rS(k)
    end do
    tmp1 = p_wrk3d*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, tmp1, rS2(1), wrk1d)
    tmp1 = p_wrk3d*tmp1
    call AVG_IK_V(imax, jmax, kmax, tmp1, rS3(1), wrk1d)
    tmp1 = p_wrk3d*tmp1
    call AVG_IK_V(imax, jmax, kmax, tmp1, rS4(1), wrk1d)

    if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        fS2(:) = rS2(:)
        fS3(:) = rS3(:)
        fS4(:) = rS4(:)

    else
        do k = 1, kmax
            p_wrk3d(:, :, k) = s_local(:, :, k) - fS(k)
        end do
        tmp1 = p_wrk3d*p_wrk3d*rho
        call AVG_IK_V(imax, jmax, kmax, tmp1, fS2(1), wrk1d)
        tmp1 = p_wrk3d*tmp1
        call AVG_IK_V(imax, jmax, kmax, tmp1, fS3(1), wrk1d)
        tmp1 = p_wrk3d*tmp1
        call AVG_IK_V(imax, jmax, kmax, tmp1, fS4(1), wrk1d)
        fS2(:) = fS2(:)/rR(:)
        fS3(:) = fS3(:)/rR(:)
        fS4(:) = fS4(:)/rR(:)

    end if

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fS2(1), Rss_z(1))

    ! -----------------------------------------------------------------------
    ! Turbulent fluxes
    do k = 1, kmax
        p_wrk3d(:, :, k) = s_local(:, :, k) - fS(k)
    end do
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) p_wrk3d = p_wrk3d*rho

    do k = 1, kmax
        dsdx(:, :, k) = p_wrk3d(:, :, k)*(u(:, :, k) - fU(k))
        dsdy(:, :, k) = p_wrk3d(:, :, k)*(v(:, :, k) - fV(k))
        dsdz(:, :, k) = p_wrk3d(:, :, k)*(w(:, :, k) - fW(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, dsdx, Rsu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Rsv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Rsw(1), wrk1d)
    Rsu(:) = Rsu(:)/rR(:)
    Rsv(:) = Rsv(:)/rR(:)
    Rsw(:) = Rsw(:)/rR(:)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rsu(1), Rsu_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rsv(1), Rsv_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rsw(1), Rsw_z(1))

    ! -----------------------------------------------------------------------
    ! Turbulent transport terms
    do k = 1, kmax
        tmp1(:, :, k) = dsdy(:, :, k)*(s_local(:, :, k) - fS(k))
        dsdx(:, :, k) = dsdx(:, :, k)*(w(:, :, k) - fW(k))
        dsdy(:, :, k) = dsdy(:, :, k)*(w(:, :, k) - fW(k))
        dsdz(:, :, k) = dsdz(:, :, k)*(w(:, :, k) - fW(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, tmp1, Tssz1(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdx, Tsuz1(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Tsvz1(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Tswz1(1), wrk1d)

    ! -----------------------------------------------------------------------
    ! Pressure terms in transport equations
    call AVG_IK_V(imax, jmax, kmax, p_loc, rP(1), wrk1d)

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s_local, dsdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s_local, dsdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s_local, dsdz)
    do k = 1, kmax
        tmp1(:, :, k) = (p_loc(:, :, k) - rP(k))*(s_local(:, :, k) - fS(k))
        dsdx(:, :, k) = (p_loc(:, :, k) - rP(k))*dsdx(:, :, k)
        dsdy(:, :, k) = (p_loc(:, :, k) - rP(k))*dsdy(:, :, k)
        dsdz(:, :, k) = (p_loc(:, :, k) - rP(k))*(dsdz(:, :, k) - fS_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, tmp1, Tswz3(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdx, PIsu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, PIsv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, PIsw(1), wrk1d)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rP(1), aux(1))
    Gsw(:) = (rS(:) - fS(:))*aux(:)

    ! #######################################################################
    ! Cross-scalar terms
    ! #######################################################################
    iv = ig(8) - 1

    do is_loc = 1, inb_scal_array
        call AVG_IK_V(imax, jmax, kmax, s(:, :, :, is_loc), aux(1), wrk1d)
        do k = 1, kmax
            tmp1(:, :, k) = (s(:, :, k, is_loc) - aux(k))*(s_local(:, :, k) - fS(k))
            tmp2(:, :, k) = tmp1(:, :, k)*(s_local(:, :, k) - fS(k))
        end do
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, tmp1, mean2d(1, iv), wrk1d)
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, tmp2, mean2d(1, iv), wrk1d)
    end do

    ! #######################################################################
    ! Source terms
    ! #######################################################################
    iv = ig(1) + offset_sources - 1
    tmp1 = 0.0_wp                               ! Accumulating sources in tmp1
    if (infraredProps%active(is)) then          ! Radiation source in tmp1
        call Radiation_Infrared_Z(infraredProps, imax, jmax, kmax, fdm_Int0, s, tmp1, tmp2, tmp3, dsdx, dsdy, dsdz)
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, tmp1)
        end if
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, tmp1, mean2d(:, iv), wrk1d) ! source
        dsdy = dsdz - dsdy
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, dsdy, mean2d(:, iv), wrk1d) ! net flux upwards
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, dsdz, mean2d(:, iv), wrk1d) ! flux up only
    end if

    if (sedimentationProps%active(is)) then     ! Sedimentation in tmp2 and dsdz
        call Microphysics_Sedimentation(sedimentationProps, imax, jmax, kmax, is, g(3), s, tmp2, tmp3, dsdx)
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, tmp2)
        end if
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, tmp2, mean2d(:, iv), wrk1d) ! source
        iv = iv + 1; call AVG_IK_V(imax, jmax, kmax, dsdx, mean2d(:, iv), wrk1d) ! net flux upwards
        tmp1 = tmp1 + tmp2
    end if

    if (evaporationProps%active(is)) then       ! Evaporation in tmp3
        ! To be done
        tmp1 = tmp1 + tmp3
    end if

    ! Total sources
    call AVG_IK_V(imax, jmax, kmax, tmp1, rQ(1), wrk1d)
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) tmp1 = tmp1*rho
    call AVG_IK_V(imax, jmax, kmax, tmp1, fQ(1), wrk1d)
    fQ(:) = fQ(:)/rR(:)

    do k = 1, kmax
        tmp1(:, :, k) = (s_local(:, :, k) - fS(k))*tmp1(:, :, k)
        dsdx(:, :, k) = (u(:, :, k) - fU(k))*tmp1(:, :, k)
        dsdy(:, :, k) = (v(:, :, k) - fV(k))*tmp1(:, :, k)
        dsdz(:, :, k) = (w(:, :, k) - fW(k))*tmp1(:, :, k)
    end do
    call AVG_IK_V(imax, jmax, kmax, tmp1, Qss(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdx, Qsu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Qsv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Qsw(1), wrk1d)
    Qss(:) = Qss(:)*2.0_wp/rR(:)
    Qsu(:) = Qsu(:)/rR(:)
    Qsv(:) = Qsv(:)/rR(:)
    Qsw(:) = Qsw(:)/rR(:)

    ! #######################################################################
    ! Derivatives
    ! #######################################################################
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), s_local, dsdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), s_local, dsdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s_local, dsdz)

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdx*dsdx + dsdy*dsdy + dsdz*dsdz
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Ess(1), wrk1d)
    Ess(:) = Ess(:)*diff*2.0_wp

    ! -----------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), u, tmp1)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), v, tmp2)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), w, tmp3)

    ! Transport term
    p_wrk3d = (tmp3*2.0_wp - tmp1 - tmp2)*c23*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_zz(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = -(p_wrk3d(:, :, k) - Tau_zz(k))*(s_local(:, :, k) - fS(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tswz2(1), wrk1d)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tau_zz(1), Tau_zz_z(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdx*((tmp1*2.0_wp - tmp2 - tmp3)*c23*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Esu(1), wrk1d)

    p_wrk3d = dsdy*((tmp2*2.0_wp - tmp1 - tmp3)*c23*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Esv(1), wrk1d)

    p_wrk3d = dsdz*((tmp3*2.0_wp - tmp1 - tmp2)*c23*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Esw(1), wrk1d)

    ! -----------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), w, tmp2)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), u, tmp1)

    ! Transport term
    p_wrk3d = (tmp1 + tmp2)*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_zx(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = -(p_wrk3d(:, :, k) - Tau_zx(k))*(s_local(:, :, k) - fS(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tsuz2(1), wrk1d)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tau_zx(1), Tau_zx_z(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdz*((tmp1 + tmp2)*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esu(:) = Esu(:) + aux(:)

    p_wrk3d = dsdx*((tmp1 + tmp2)*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esw(:) = Esw(:) + aux(:)

    ! -----------------------------------------------------------------------
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), w, tmp3)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), v, tmp2)

    ! Transport term
    p_wrk3d = (tmp3 + tmp2)*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yz(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = -(p_wrk3d(:, :, k) - Tau_yz(k))*(s_local(:, :, k) - fS(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tsvz2(1), wrk1d)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tau_yz(1), Tau_yz_z(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdz*((tmp3 + tmp2)*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esv(:) = Esv(:) + aux(:)

    p_wrk3d = dsdy*((tmp3 + tmp2)*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esw(:) = Esw(:) + aux(:)

    ! -----------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), v, tmp3)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), u, tmp1)

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdy*((tmp3 + tmp1)*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esu(:) = Esu(:) + aux(:)

    p_wrk3d = dsdx*((tmp3 + tmp1)*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esv(:) = Esv(:) + aux(:)

    ! -----------------------------------------------------------------------
    ! Moments
    do k = 1, kmax
        p_wrk3d(:, :, k) = dsdz(:, :, k) - rS_z(k)
    end do

    tmp1 = dsdx*dsdx
    tmp2 = dsdy*dsdy
    tmp3 = p_wrk3d*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, tmp1, S_x2(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp2, S_y2(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp3, S_z2(1), wrk1d)

    tmp1 = tmp1*dsdx
    tmp2 = tmp2*dsdy
    tmp3 = tmp3*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, tmp1, S_x3(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp2, S_y3(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp3, S_z3(1), wrk1d)

    tmp1 = tmp1*dsdx
    tmp2 = tmp2*dsdy
    tmp3 = tmp3*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, tmp1, S_x4(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp2, S_y4(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp3, S_z4(1), wrk1d)

    ! -----------------------------------------------------------------------
    ! Molecular fluxes
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) dsdz = dsdz*vis
    call AVG_IK_V(imax, jmax, kmax, dsdz, Fz(1), wrk1d)

    ! Contribution to turbulent transport
    do k = 1, kmax
        p_wrk3d(:, :, k) = (dsdz(:, :, k) - Fz(k))*(s_local(:, :, k) - fS(k))
        tmp1(:, :, k) = (dsdz(:, :, k) - Fz(k))*(u(:, :, k) - fU(k))
        tmp2(:, :, k) = (dsdz(:, :, k) - Fz(k))*(v(:, :, k) - fV(k))
        tmp3(:, :, k) = (dsdz(:, :, k) - Fz(k))*(w(:, :, k) - fW(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tssz2(1), wrk1d)
    Tssz2(:) = -Tssz2(:)*diff*2.0_wp
    call AVG_IK_V(imax, jmax, kmax, tmp1, aux(1), wrk1d)
    Tsuz2(:) = Tsuz2(:) - aux(:)*diff
    call AVG_IK_V(imax, jmax, kmax, tmp2, aux(1), wrk1d)
    Tsvz2(:) = Tsvz2(:) - aux(:)*diff
    call AVG_IK_V(imax, jmax, kmax, tmp3, aux(1), wrk1d)
    Tswz2(:) = Tswz2(:) - aux(:)*diff

    Fz(:) = Fz(:)*diff
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Fz(1), Fz_z(1))

    ! Contribution to dissipation
    Ess(:) = (Ess(:) - Fz(:)*rS_z(:) - Fz(:)*rS_z(:))/rR(:)
    Esu(:) = (Esu(:) - Tau_zx(:)*rS_z(:) - Fz(:)*rU_z(:))/rR(:)
    Esv(:) = (Esv(:) - Tau_yz(:)*rS_z(:) - Fz(:)*rV_z(:))/rR(:)
    Esw(:) = (Esw(:) - Tau_zz(:)*rS_z(:) - Fz(:)*rW_z(:))/rR(:)

    ! #######################################################################
    ! Source terms in transport equations
    ! #######################################################################
    select case (nse_eqns)
    case (DNS_EQNS_ANELASTIC)
        ! call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, p_wrk3d)

    case (DNS_EQNS_BOUSSINESQ)
        call Gravity_Source(gravityProps, imax, jmax, kmax, s, p_wrk3d)

    case (DNS_EQNS_COMPRESSIBLE)
        p_wrk3d = rho

    end select
    do k = 1, kmax
        p_wrk3d(:, :, k) = (s_local(:, :, k) - fS(k))*p_wrk3d(:, :, k)
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Bsu(:) = aux(:)*gravityProps%vector(1)/rR(:)
    Bsv(:) = aux(:)*gravityProps%vector(2)/rR(:)
    Bsw(:) = aux(:)*gravityProps%vector(3)/rR(:)

    ! #######################################################################
    ! Complete budget equations
    ! #######################################################################
    ! Transport terms
    aux(:) = Tssz1(:) + Tssz2(:)
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), aux(1), Tssz_z(1))
    aux(:) = Tsuz1(:) + Tsuz2(:)
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), aux(1), Tsuz_z(1))
    aux(:) = Tsvz1(:) + Tsvz2(:)
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), aux(1), Tsvz_z(1))
    aux(:) = Tswz1(:) + Tswz2(:) + Tswz3(:)
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), aux(1), Tswz_z(1))

    ! Convective terms
    Css(:) = -fW(:)*Rss_z(:)
    Csu(:) = -fW(:)*Rsu_z(:)
    Csv(:) = -fW(:)*Rsv_z(:)
    Csw(:) = -fW(:)*Rsw_z(:)

    ! Production terms
    Pss(:) = -Rsw(:)*fS_z(:)*2.0_wp
    Psu(:) = -Rsw(:)*fU_z(:) - Rwu(:)*fS_z(:)
    Psv(:) = -Rsw(:)*fV_z(:) - Rwv(:)*fS_z(:)
    Psw(:) = -Rsw(:)*fW_z(:) - Rww(:)*fS_z(:)

    ! Diffusion variable-density terms
    Dss(:) = (rS(:) - fS(:))*Fz_z(:)*2.0_wp
    Dsu(:) = (rS(:) - fS(:))*Tau_zx_z(:) + (rU(:) - fU(:))*Fz_z(:)
    Dsv(:) = (rS(:) - fS(:))*Tau_yz_z(:) + (rV(:) - fV(:))*Fz_z(:)
    Dsw(:) = (rS(:) - fS(:))*Tau_zz_z(:) + (rW(:) - fW(:))*Fz_z(:)

    ! ! Coriolis terms
    ! dummy = coriolis%vector(2)
    ! Fsu(:) = dummy*Rsw(:)
    ! Fsw(:) = -dummy*Rsu(:)

    ! Transient terms
    Rss_t(:) = Css(:) + Pss(:) - Ess(:) + Qss(:) + (Dss(:) - Tssz_z(:))/rR(:)
    Rsu_t(:) = Csu(:) + Psu(:) - Esu(:) + Bsu(:) - Fsu(:) + Qsu(:) + (PIsu(:) + Dsu(:) - Gsu(:) - Tsuz_z(:))/rR(:)
    Rsv_t(:) = Csv(:) + Psv(:) - Esv(:) + Bsv(:) - Fsv(:) + Qsv(:) + (PIsv(:) + Dsv(:) - Gsv(:) - Tsvz_z(:))/rR(:)
    Rsw_t(:) = Csw(:) + Psw(:) - Esw(:) + Bsw(:) - Fsw(:) + Qsw(:) + (PIsw(:) + Dsw(:) - Gsw(:) - Tswz_z(:))/rR(:)

    ! ###################################################################
    ! Output
    ! #######################################################################
    ! 11 t-dependent variables, for consistency with old format
    ! ng = ng +1
    ! groupname(ng) = ''
    ! varname(ng)   = 'dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy'
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''

    write (line1, *) is; line1 = 'avg'//trim(adjustl(line1))//'s'
    write (name, *) itime; name = trim(adjustl(line1))//trim(adjustl(name))
    call IO_WRITE_AVERAGES(name, itime, rtime, kmax, nv, ng, z%nodes, varname, groupname, mean2d)

    return
end subroutine AVG_SCAL_XZ
