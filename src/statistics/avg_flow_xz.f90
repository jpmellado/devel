#include "tlab_error.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOy, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array dudz until pressure block
!#
!# Reynolds and Favre averages
!#
!########################################################################

subroutine AVG_FLOW_XZ(q, s, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, mean2d)
    use TLab_Constants, only: wp, wi, MAX_AVG_TEMPORAL
    use TLab_Constants, only: efile, lfile
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
#endif
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_flow_array, inb_scal_array
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: u, v, w, p_wrk3d
    use Thermo_Compressible, only: rho, vis
    use TLab_Grid, only: z
    use FDM, only: g
    use OPR_Partial
    use Thermodynamics, only: imode_thermo, THERMO_TYPE_NONE, THERMO_TYPE_ANELASTIC, THERMO_TYPE_COMPRESSIBLE
    use Thermo_Base, only: Thermo_Psat_Polynomial
    use Thermo_AirWater, only: inb_scal_T
    use Thermo_Anelastic
    use NavierStokes
    use Gravity, only: gravityProps, Gravity_Source
    ! use Rotation, only: coriolis
    use Averages, only: AVG_IK_V

    implicit none

    real(wp), intent(in) :: q(imax, jmax, kmax, inb_flow_array)
    real(wp), intent(in) :: s(imax, jmax, kmax, inb_scal_array)
    real(wp), intent(inout) :: dudx(imax, jmax, kmax)
    real(wp), intent(inout) :: dudy(imax, jmax, kmax)
    real(wp), intent(inout) :: dudz(imax, jmax, kmax)
    real(wp), intent(inout) :: dvdx(imax, jmax, kmax)
    real(wp), intent(inout) :: dvdy(imax, jmax, kmax)
    real(wp), intent(inout) :: dvdz(imax, jmax, kmax)
    real(wp), intent(inout) :: dwdx(imax, jmax, kmax)
    real(wp), intent(inout) :: dwdy(imax, jmax, kmax)
    real(wp), intent(inout) :: dwdz(imax, jmax, kmax)
    real(wp), intent(inout) :: mean2d(kmax, MAX_AVG_TEMPORAL)

    target q, dudz

    ! -------------------------------------------------------------------
    integer, parameter :: MAX_VARS_GROUPS = 20
    integer k
    real(wp) dummy
    real(wp) c23

    integer ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nv

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
        p_loc => dudz
    end if

    c23 = 2.0_wp/3.0_wp

    ! Variable definition and memory management
    ! -----------------------------------------------------------------------
    ng = 1; ig(ng) = 1
#define rU(k)     mean2d(k,ig(1))
#define rV(k)     mean2d(k,ig(1)+1)
#define rW(k)     mean2d(k,ig(1)+2)
#define fU(k)     mean2d(k,ig(1)+3)
#define fV(k)     mean2d(k,ig(1)+4)
#define fW(k)     mean2d(k,ig(1)+5)
#define rP(k)     mean2d(k,ig(1)+6)

#define Rxx(k)    mean2d(k,ig(1)+7)
#define Ryy(k)    mean2d(k,ig(1)+8)
#define Rzz(k)    mean2d(k,ig(1)+9)
#define Rxy(k)    mean2d(k,ig(1)+10)
#define Rxz(k)    mean2d(k,ig(1)+11)
#define Ryz(k)    mean2d(k,ig(1)+12)
#define Tke(k)    mean2d(k,ig(1)+13)
#define rP2(k)    mean2d(k,ig(1)+14)
    sg(ng) = 15

    groupname(ng) = 'Dynamics'
    varname(ng) = 'rU rV rW fU fV fW rP Rxx Ryy Rzz Rxy Rxz Ryz Tke rP2'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define vortx(k)  mean2d(k,ig(2)  )
#define vorty(k)  mean2d(k,ig(2)+1)
#define vortz(k)  mean2d(k,ig(2)+2)
#define vortx2(k) mean2d(k,ig(2)+3)
#define vorty2(k) mean2d(k,ig(2)+4)
#define vortz2(k) mean2d(k,ig(2)+5)
    sg(ng) = 6

    groupname(ng) = 'Vorticity'
    varname(ng) = 'Wx Wy Wz Wx2 Wy2 Wz2'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxx_t(k)  mean2d(k,ig(3)  )
#define Bxx(k)    mean2d(k,ig(3)+1)
#define Cxx(k)    mean2d(k,ig(3)+2)
#define Pxx(k)    mean2d(k,ig(3)+3)
#define Exx(k)    mean2d(k,ig(3)+4)
#define PIxx(k)   mean2d(k,ig(3)+5)
#define Fxx(k)    mean2d(k,ig(3)+6)
#define Txxz_z(k) mean2d(k,ig(3)+7)
#define Txxz(k)   mean2d(k,ig(3)+8)
#define Gxx(k)    mean2d(k,ig(3)+9)
#define Dxx(k)    mean2d(k,ig(3)+10)
    sg(ng) = 11

    groupname(ng) = 'RxxBudget'
    varname(ng) = 'Rxx_t Bxx Cxx Pxx Exx PIxx Fxx Txxz_z Txxz Gxx Dxx'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Ryy_t(k)  mean2d(k,ig(4)  )
#define Byy(k)    mean2d(k,ig(4)+1)
#define Cyy(k)    mean2d(k,ig(4)+2)
#define Pyy(k)    mean2d(k,ig(4)+3)
#define Eyy(k)    mean2d(k,ig(4)+4)
#define PIyy(k)   mean2d(k,ig(4)+5)
#define Fyy(k)    mean2d(k,ig(4)+6)
#define Tyyz_z(k) mean2d(k,ig(4)+7)
#define Tyyz(k)   mean2d(k,ig(4)+8)
#define Gyy(k)    mean2d(k,ig(4)+9)
#define Dyy(k)    mean2d(k,ig(4)+10)
    sg(ng) = 11

    groupname(ng) = 'RyyBudget'
    varname(ng) = 'Ryy_t Byy Cyy Pyy Eyy PIyy Fyy Tyyz_z Tyyz Gyy Dyy'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rzz_t(k)  mean2d(k,ig(5)  )
#define Bzz(k)    mean2d(k,ig(5)+1)
#define Czz(k)    mean2d(k,ig(5)+2)
#define Pzz(k)    mean2d(k,ig(5)+3)
#define Ezz(k)    mean2d(k,ig(5)+4)
#define PIzz(k)   mean2d(k,ig(5)+5)
#define Fzz(k)    mean2d(k,ig(5)+6)
#define Tzzz_z(k) mean2d(k,ig(5)+7)
#define Tzzz(k)   mean2d(k,ig(5)+8)
#define Gzz(k)    mean2d(k,ig(5)+9)
#define Dzz(k)    mean2d(k,ig(5)+10)
    sg(ng) = 11

    groupname(ng) = 'RzzBudget'
    varname(ng) = 'Rzz_t Bzz Czz Pzz Ezz PIzz Fzz Tzzz_z Tzzz Gzz Dzz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxy_t(k)  mean2d(k,ig(6)  )
#define Bxy(k)    mean2d(k,ig(6)+1)
#define Cxy(k)    mean2d(k,ig(6)+2)
#define Pxy(k)    mean2d(k,ig(6)+3)
#define Exy(k)    mean2d(k,ig(6)+4)
#define PIxy(k)   mean2d(k,ig(6)+5)
#define Fxy(k)    mean2d(k,ig(6)+6)
#define Txyz_z(k) mean2d(k,ig(6)+7)
#define Txyz(k)   mean2d(k,ig(6)+8)
#define Gxy(k)    mean2d(k,ig(6)+9)
#define Dxy(k)    mean2d(k,ig(6)+10)
    sg(ng) = 11

    groupname(ng) = 'RxyBudget'
    varname(ng) = 'Rxy_t Bxy Cxy Pxy Exy PIxy Fxy Txyz_z Txyz Gxy Dxy'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxz_t(k)  mean2d(k,ig(7)  )
#define Bxz(k)    mean2d(k,ig(7)+1)
#define Cxz(k)    mean2d(k,ig(7)+2)
#define Pxz(k)    mean2d(k,ig(7)+3)
#define Exz(k)    mean2d(k,ig(7)+4)
#define PIxz(k)   mean2d(k,ig(7)+5)
#define Fxz(k)    mean2d(k,ig(7)+6)
#define Txzz_z(k) mean2d(k,ig(7)+7)
#define Txzz(k)   mean2d(k,ig(7)+8)
#define Gxz(k)    mean2d(k,ig(7)+9)
#define Dxz(k)    mean2d(k,ig(7)+10)
    sg(ng) = 11

    groupname(ng) = 'RxzBudget'
    varname(ng) = 'Rxz_t Bxz Cxz Pxz Exz PIxz Fxz Txzz_z Txzz Gxz Dxz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Ryz_t(k)  mean2d(k,ig(8)  )
#define Byz(k)    mean2d(k,ig(8)+1)
#define Cyz(k)    mean2d(k,ig(8)+2)
#define Pyz(k)    mean2d(k,ig(8)+3)
#define Eyz(k)    mean2d(k,ig(8)+4)
#define PIyz(k)   mean2d(k,ig(8)+5)
#define Fyz(k)    mean2d(k,ig(8)+6)
#define Tyzz_z(k) mean2d(k,ig(8)+7)
#define Tyzz(k)   mean2d(k,ig(8)+8)
#define Gyz(k)    mean2d(k,ig(8)+9)
#define Dyz(k)    mean2d(k,ig(8)+10)
    sg(ng) = 11

    groupname(ng) = 'RyzBudget'
    varname(ng) = 'Ryz_t Byz Cyz Pyz Eyz PIyz Fyz Tyzz_z Tyzz Gyz Dyz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Tke_t(k)  mean2d(k,ig(9)  )
#define Buo(k)    mean2d(k,ig(9)+1)
#define Con(k)    mean2d(k,ig(9)+2)
#define Prd(k)    mean2d(k,ig(9)+3)
#define Eps(k)    mean2d(k,ig(9)+4)
#define Pi(k)     mean2d(k,ig(9)+5)
#define Tz_z(k)   mean2d(k,ig(9)+6)
#define Tz1(k)    mean2d(k,ig(9)+7)
#define Tz2(k)    mean2d(k,ig(9)+8)
#define Tz3(k)    mean2d(k,ig(9)+9)
#define Tz1_z(k)  mean2d(k,ig(9)+10)
#define Tz2_z(k)  mean2d(k,ig(9)+11)
#define Tz3_z(k)  mean2d(k,ig(9)+12)
#define Gkin(k)   mean2d(k,ig(9)+13)
#define Dkin(k)   mean2d(k,ig(9)+14)
#define Phi(k)    mean2d(k,ig(9)+15)
#define ugradp(k) mean2d(k,ig(9)+16)
    sg(ng) = 17

    groupname(ng) = 'TkeBudget'
    varname(ng) = 'Tke_t Buo Con Prd Eps Pi Trp Trp1 Trp2 Trp3 Trp1_z Trp2_z Trp3_z G D Phi UgradP'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rU3(k)   mean2d(k,ig(10)  )
#define rU4(k)   mean2d(k,ig(10)+1)
#define rV3(k)   mean2d(k,ig(10)+2)
#define rV4(k)   mean2d(k,ig(10)+3)
#define rW3(k)   mean2d(k,ig(10)+4)
#define rW4(k)   mean2d(k,ig(10)+5)
    sg(ng) = 6

    groupname(ng) = 'HigherOrder'
    varname(ng) = 'rU3 rU4 rV3 rV4 rW3 rW4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define U_z1(k)  mean2d(k,ig(11)  )
#define V_z1(k)  mean2d(k,ig(11)+1)
#define W_z1(k)  mean2d(k,ig(11)+2)
#define U_ii2(k) mean2d(k,ig(11)+3)
#define U_x2(k)  mean2d(k,ig(11)+4)
#define U_y2(k)  mean2d(k,ig(11)+5)
#define U_z2(k)  mean2d(k,ig(11)+6)
#define V_x2(k)  mean2d(k,ig(11)+7)
#define V_y2(k)  mean2d(k,ig(11)+8)
#define V_z2(k)  mean2d(k,ig(11)+9)
#define W_x2(k)  mean2d(k,ig(11)+10)
#define W_y2(k)  mean2d(k,ig(11)+11)
#define W_z2(k)  mean2d(k,ig(11)+12)
#define U_x3(k)  mean2d(k,ig(11)+13)
#define U_y3(k)  mean2d(k,ig(11)+14)
#define U_z3(k)  mean2d(k,ig(11)+15)
#define V_x3(k)  mean2d(k,ig(11)+16)
#define V_y3(k)  mean2d(k,ig(11)+17)
#define V_z3(k)  mean2d(k,ig(11)+18)
#define W_x3(k)  mean2d(k,ig(11)+19)
#define W_y3(k)  mean2d(k,ig(11)+20)
#define W_z3(k)  mean2d(k,ig(11)+21)
#define U_x4(k)  mean2d(k,ig(11)+22)
#define U_y4(k)  mean2d(k,ig(11)+23)
#define U_z4(k)  mean2d(k,ig(11)+24)
#define V_x4(k)  mean2d(k,ig(11)+25)
#define V_y4(k)  mean2d(k,ig(11)+26)
#define V_z4(k)  mean2d(k,ig(11)+27)
#define W_x4(k)  mean2d(k,ig(11)+28)
#define W_y4(k)  mean2d(k,ig(11)+29)
#define W_z4(k)  mean2d(k,ig(11)+30)
    sg(ng) = 31

    groupname(ng) = 'DerivativeFluctuations'
    varname(ng) = 'U_z1 V_z1 W_z1 U_ii2 ' &
                  //'U_x2 U_y2 U_z2 V_x2 V_y2 V_z2 W_x2 W_y2 W_z2 ' &
                  //'U_x3 U_y3 U_z3 V_x3 V_y3 V_z3 W_x3 W_y3 W_z3 ' &
                  //'U_x4 U_y4 U_z4 V_x4 V_y4 V_z4 W_x4 W_y4 W_z4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rB(k)        mean2d(k,ig(12)  )
#define Pot(k)       mean2d(k,ig(12)+1)
#define bfreq_fr(k)  mean2d(k,ig(12)+2)
#define bfreq_eq(k)  mean2d(k,ig(12)+3)
    sg(ng) = 4

    groupname(ng) = 'Buoyancy'
    varname(ng) = 'rB Pot BuoyFreq_fr BuoyFreq_eq'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define pref(k)      mean2d(k,ig(13))
#define rref(k)      mean2d(k,ig(13)+1)
#define tref(k)      mean2d(k,ig(13)+2)
#define rR(k)        mean2d(k,ig(13)+3)
#define rT(k)        mean2d(k,ig(13)+4)
#define rR2(k)       mean2d(k,ig(13)+5)
#define rT2(k)       mean2d(k,ig(13)+6)
#define theta(k)     mean2d(k,ig(13)+7)
#define thetav(k)    mean2d(k,ig(13)+8)
#define thetae(k)    mean2d(k,ig(13)+9)
#define thetal(k)    mean2d(k,ig(13)+10)
#define mse(k)       mean2d(k,ig(13)+11)
#define psat(k)      mean2d(k,ig(13)+12)
#define rh(k)        mean2d(k,ig(13)+13)
#define dewpoint(k)  mean2d(k,ig(13)+14)
#define lapse_fr(k)  mean2d(k,ig(13)+15)
#define lapse_eq(k)  mean2d(k,ig(13)+16)
#define lapse_dew(k) mean2d(k,ig(13)+17)
    sg(ng) = 18

    groupname(ng) = 'Thermodynamics'
    varname(ng) = 'pRef rRef TRef rR rT rR2 rT2 ' &
                  //'Theta ThetaV ThetaE ThetaL MoistStaticEnergy ' &
                  //'Psat RelativeHumidity DewPoint ' &
                  //'LapseRate_fr LapseRate_eq LapseRate_dew'

    ! -----------------------------------------------------------------------
    ! Auxiliary variables depending on z and t; this last group is not written
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rUf(k)    mean2d(k,ig(14))
#define rVf(k)    mean2d(k,ig(14)+1)
#define rWf(k)    mean2d(k,ig(14)+2)

#define rU_z(k)   mean2d(k,ig(14)+3)
#define rV_z(k)   mean2d(k,ig(14)+4)
#define rW_z(k)   mean2d(k,ig(14)+5)
#define fU_z(k)   mean2d(k,ig(14)+6)
#define fV_z(k)   mean2d(k,ig(14)+7)
#define fW_z(k)   mean2d(k,ig(14)+8)
#define rP_z(k)   mean2d(k,ig(14)+9)

#define Rxx_z(k)  mean2d(k,ig(14)+10)
#define Ryy_z(k)  mean2d(k,ig(14)+11)
#define Rzz_z(k)  mean2d(k,ig(14)+12)
#define Rxy_z(k)  mean2d(k,ig(14)+13)
#define Rxz_z(k)  mean2d(k,ig(14)+14)
#define Ryz_z(k)  mean2d(k,ig(14)+15)

#define Tau_xz(k)   mean2d(k,ig(14)+16)
#define Tau_yz(k)   mean2d(k,ig(14)+17)
#define Tau_zz(k)   mean2d(k,ig(14)+18)
#define Tau_xz_z(k) mean2d(k,ig(14)+19)
#define Tau_yz_z(k) mean2d(k,ig(14)+20)
#define Tau_zz_z(k) mean2d(k,ig(14)+21)

#define aux(k)      mean2d(k,ig(14)+22)

    sg(ng) = 23

    ! -----------------------------------------------------------------------
    nv = ig(ng) + sg(ng) - 1
    if (MAX_AVG_TEMPORAL < nv) then
        call TLab_Write_ASCII(efile, __FILE__//'. Not enough space in local arrays.')
        call TLab_Stop(DNS_ERROR_AVGTMP)
    end if
    mean2d(:, 1:nv) = 0.0_wp

    ng = ng - 1
    nv = ig(ng) + sg(ng) - 1        ! the last group with auxiliary variables is not written out

    select case (imode_thermo)
    case (THERMO_TYPE_NONE)
        ng = ng - 1
        nv = ig(ng) + sg(ng) - 1    ! the thermodynamics group is not written out

    case (THERMO_TYPE_ANELASTIC)

    case (THERMO_TYPE_COMPRESSIBLE)

    end select

    ! #######################################################################
    write (line1, *) itime; line1 = 'Calculating flow statistics at It'//trim(adjustl(line1))//'...'
    call TLab_Write_ASCII(lfile, line1)

    ! ###################################################################
    ! Averages (do not overwrite dudz; it contains p for incompressible case)
    ! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, __FILE__//': Section 2')
#endif

    ! Velocity
    call AVG_IK_V(imax, jmax, kmax, u, rU(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, v, rV(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, w, rW(1), wrk1d)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rU(1), rU_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rV(1), rV_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rW(1), rW_z(1))

    U_z1(:) = rU_z(:)
    V_z1(:) = rV_z(:)
    W_z1(:) = rW_z(:)

    ! Density and Favre avrages
    if (nse_eqns == DNS_EQNS_BOUSSINESQ) then
        fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

    else if (nse_eqns == DNS_EQNS_ANELASTIC) then
        fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

    else
        call AVG_IK_V(imax, jmax, kmax, rho, rR(1), wrk1d)

        dwdx = rho*u
        dwdy = rho*v
        dwdz = rho*w
        call AVG_IK_V(imax, jmax, kmax, dwdx, fU(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dwdy, fV(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dwdz, fW(1), wrk1d)
        fU(:) = fU(:)/rR(:)
        fV(:) = fV(:)/rR(:)
        fW(:) = fW(:)/rR(:)

    end if

    rUf(:) = rU(:) - fU(:)
    rVf(:) = rV(:) - fV(:)
    rWf(:) = rW(:) - fW(:)

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fU(1), fU_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fV(1), fV_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), fW(1), fW_z(1))

    ! Pressure
    call AVG_IK_V(imax, jmax, kmax, p_loc, rP(1), wrk1d)
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rP(1), rP_z(1))

    ! #######################################################################
    ! Main covariances (do not overwrite dudz; it contains p for incompressible case)
    ! #######################################################################
    do k = 1, kmax
        dwdx(:, :, k) = u(:, :, k) - fU(k)
        dwdy(:, :, k) = v(:, :, k) - fV(k)
        dwdz(:, :, k) = w(:, :, k) - fW(k)
    end do

    dvdx = dwdx*dwdx
    dvdy = dwdy*dwdy
    dvdz = dwdz*dwdz
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Rxx(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Ryy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Rzz(1), wrk1d)
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        Rxx(:) = Rxx(:)/rR(:)
        Ryy(:) = Ryy(:)/rR(:)
        Rzz(:) = Rzz(:)/rR(:)
    end if

    dvdx = dwdx*dwdy
    dvdy = dwdx*dwdz
    dvdz = dwdy*dwdz
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Rxy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Rxz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Ryz(1), wrk1d)
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        Rxy(:) = Rxy(:)/rR(:)
        Rxz(:) = Rxz(:)/rR(:)
        Ryz(:) = Ryz(:)/rR(:)
    end if

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rxx(1), Rxx_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Ryy(1), Ryy_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rzz(1), Rzz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rxy(1), Rxy_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Rxz(1), Rxz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Ryz(1), Ryz_z(1))

    ! higher-order moments
    p_wrk3d = dwdx*dwdx*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rU3(1), wrk1d)
    p_wrk3d = dwdx*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rU4(1), wrk1d)

    p_wrk3d = dwdy*dwdy*dwdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rV3(1), wrk1d)
    p_wrk3d = dwdy*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rV4(1), wrk1d)

    p_wrk3d = dwdz*dwdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rW3(1), wrk1d)
    p_wrk3d = dwdz*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rW4(1), wrk1d)

    ! Triple-velocity correlations
    dvdx = dwdx*dwdx*dwdz
    dvdy = dwdy*dwdy*dwdz
    dvdz = dwdz*dwdz*dwdz
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Txxz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Tyyz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Tzzz(1), wrk1d)
    Tz1(:) = (Txxz(:) + Tyyz(:) + Tzzz(:))*0.5_wp

    dvdx = dwdx*dwdy*dwdz
    dvdy = dwdx*dwdz*dwdz
    dvdz = dwdy*dwdz*dwdz
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Txyz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Txzz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Tyzz(1), wrk1d)

    ! Pressure
    do k = 1, kmax
        dvdz(:, :, k) = p_loc(:, :, k) - rP(k)
    end do
    p_wrk3d = dvdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rP2(1), wrk1d)

    ! Pressure-velocity correlation in TKE transport terms
    dwdx = dwdx*dvdz
    dwdy = dwdy*dvdz
    dwdz = dwdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, dwdx, aux(1), wrk1d)
    Txzz(:) = Txzz(:) + aux(:)
    call AVG_IK_V(imax, jmax, kmax, dwdy, aux(1), wrk1d)
    Tyzz(:) = Tyzz(:) + aux(:)
    call AVG_IK_V(imax, jmax, kmax, dwdz, Tz2(1), wrk1d)
    Tzzz(:) = Tzzz(:) + Tz2(:)*2.0_wp

    ! ###################################################################
    ! Pressure; array dudz containing p is used only up to this section
    !
    ! dudx = du/dx
    ! dudy = du/dy
    ! dudz = p
    ! dvdx = dv/dx
    ! dvdy = dv/dy
    ! dvdz = p_prime
    ! dwdx =       ; dp/dx
    ! dwdy =       ; dp/dy
    ! dwdz = dw/dz ; dp/dz
    ! ###################################################################
    ! Pressure convection term
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), p_loc, dwdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), p_loc, dwdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), p_loc, dwdz)
    p_wrk3d = u*dwdx + v*dwdy + w*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, ugradp(1), wrk1d)

    ! Pressure Strain Terms
    ! 9 derivatives are here recomputed; ok, this routine is not called that often
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), u, dudx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), v, dvdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), w, dwdz)
    dudx = dvdz*dudx ! dvdz contains the pressure fluctuation
    dvdy = dvdz*dvdy
    dwdz = dvdz*dwdz ! no need to substract rW_z
    call AVG_IK_V(imax, jmax, kmax, dudx, PIxx(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, PIyy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dwdz, PIzz(1), wrk1d)
    PIxx(:) = PIxx(:)*2.0_wp
    PIyy(:) = PIyy(:)*2.0_wp
    PIzz(:) = PIzz(:)*2.0_wp

    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), u, dudy)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), v, dvdx)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), u, dwdz) !dudz not free
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), w, dwdx)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), v, dudx) !dvdz not free
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), w, dvdy)
    dudy = dvdz*(dudy + dvdx)
    dwdz = dvdz*(dwdz + dwdx) ! no need to substract rU_z
    dudx = dvdz*(dudx + dvdy) ! no need to substract rW_z
    call AVG_IK_V(imax, jmax, kmax, dudy, PIxy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dwdz, PIxz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dudx, PIyz(1), wrk1d)

    ! ###################################################################
    ! Thermodynamic variables
    ! ###################################################################
    select case (imode_thermo)
    case (THERMO_TYPE_NONE)
    case (THERMO_TYPE_ANELASTIC)
        pref(:) = pbackground(:)
        rref(:) = rbackground(:)
        tref(:) = tbackground(:)

        call Thermo_Anelastic_Rho(imax, jmax, kmax, s, dwdx, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dwdx, rR(1), wrk1d)
        do k = 1, kmax
            dvdx(:, :, k) = (dwdx(:, :, k) - rR(k))**2
        end do
        call AVG_IK_V(imax, jmax, kmax, dvdx, rR2(1), wrk1d)

        call AVG_IK_V(imax, jmax, kmax, s(:, :, :, inb_scal_T), rT(1), wrk1d)
        do k = 1, kmax
            dvdx(:, :, k) = (s(:, :, k, inb_scal_T) - rT(k))**2
        end do
        call AVG_IK_V(imax, jmax, kmax, dvdx, rT2(1), wrk1d)

        call Thermo_Anelastic_Theta(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, theta(1), wrk1d)

        call Thermo_Anelastic_ThetaV(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, thetav(1), wrk1d)

        call Thermo_Anelastic_ThetaL(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, thetal(1), wrk1d)

        call Thermo_Anelastic_ThetaE(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, thetae(1), wrk1d)

        call Thermo_Anelastic_MSE(imax, jmax, kmax, s, dvdz)
        call AVG_IK_V(imax, jmax, kmax, dvdz, mse(1), wrk1d)

        call Thermo_Psat_Polynomial(imax*jmax*kmax, s(:, :, :, inb_scal_T), dvdz)
        call AVG_IK_V(imax, jmax, kmax, dvdz, psat(1), wrk1d)

        call Thermo_Anelastic_RH(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, rh(1), wrk1d)

        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), s(:, :, :, inb_scal_T), dudz)

        call Thermo_Anelastic_LapseEquilibrium(imax, jmax, kmax, s, dwdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dwdz, lapse_eq(1), wrk1d)
        !frequency(:) = (lapse(:) + dTdy(:))/locT(:)
        p_wrk3d = (dwdz + dudz)/dwdx
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, bfreq_eq(1), wrk1d)
        bfreq_eq(:) = bfreq_eq(:)*gravityProps%vector(3)

        call Thermo_Anelastic_LapseFrozen(imax, jmax, kmax, s, dwdz)
        call AVG_IK_V(imax, jmax, kmax, dwdz, lapse_fr(1), wrk1d)
        !frequency(:) = (lapse(:) + dTdy(:))/locT(:)
        p_wrk3d = (dwdz + dudz)/dwdx
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, bfreq_fr(1), wrk1d)
        bfreq_fr(:) = bfreq_fr(:)*gravityProps%vector(3)

        call Thermo_Anelastic_Pvapor(imax, jmax, kmax, s, dudz)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), dudz, dudy)
        ! dwdz should contains lapse_fr, since lapse_dew = lapse_fr when saturated
        call Thermo_Anelastic_Weight_DewPoint(imax, jmax, kmax, s, dudy, p_wrk3d, dwdz)
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, dewpoint(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dwdz, lapse_dew(1), wrk1d)

    end select

    ! ###################################################################
    ! Buoyancy
    ! ###################################################################
    if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then

        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ)
            call Gravity_Source(gravityProps, imax, jmax, kmax, s, dudx)

            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), dudx, dudy)
            call AVG_IK_V(imax, jmax, kmax, dudy, bfreq_fr(1), wrk1d)
            bfreq_fr(:) = bfreq_fr(:)*gravityProps%vector(3)
            bfreq_eq(:) = bfreq_fr(:)

        case (DNS_EQNS_ANELASTIC)
            call Thermo_Anelastic_Buoyancy(imax, jmax, kmax, s, dudx)

        end select

        call AVG_IK_V(imax, jmax, kmax, dudx, rB(1), wrk1d)
        do k = 1, kmax
            dvdx(:, :, k) = (u(:, :, k) - rU(k))*(dudx(:, :, k) - rB(k))
            dvdy(:, :, k) = (v(:, :, k) - rV(k))*(dudx(:, :, k) - rB(k))
            dvdz(:, :, k) = (w(:, :, k) - rW(k))*(dudx(:, :, k) - rB(k))
        end do
        call AVG_IK_V(imax, jmax, kmax, dvdx, Bxx(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdy, Byy(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, Bzz(1), wrk1d)

        dummy = 1.0_wp/froude
        rB(:) = rB(:)*dummy
        ! call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), rB(1), rB_z(1))

        ! pmod(:) = -rP_z(:) + sign(rB(:), gravityProps%vector(2))

        ! potential energy
        Pot(:) = -rB(:)*z%nodes(:)

    else
        Bxx(:) = -rR(:)*rUf(:)
        Byy(:) = -rR(:)*rVf(:)
        Bzz(:) = -rR(:)*rWf(:)

        ! pmod(:) = -rP_z(:) + gravityProps%vector(2)*rR(:)

        ! potential energy
        Pot(:) = -rR(:)*z%nodes(:)*gravityProps%vector(3)

    end if

    ! gravityProps%vector includes the Froude
    Bxy(:) = Bxx(:)*gravityProps%vector(2) + Byy(:)*gravityProps%vector(1)
    Bxz(:) = Bxx(:)*gravityProps%vector(3) + Bzz(:)*gravityProps%vector(1)
    Byz(:) = Byy(:)*gravityProps%vector(3) + Bzz(:)*gravityProps%vector(2)

    Bxx(:) = 2.0_wp*Bxx(:)*gravityProps%vector(1)
    Byy(:) = 2.0_wp*Byy(:)*gravityProps%vector(2)
    Bzz(:) = 2.0_wp*Bzz(:)*gravityProps%vector(3)

    ! ###################################################################
    ! # Array storage of velocity gradient tensor
    ! #
    ! # dudx = d U / d x
    ! # dudy = d U / d y
    ! # dudz = d U / d z
    ! # dvdx = d V / d x
    ! # dvdy = d V / d y
    ! # dvdz = d V / d z
    ! # dwdx = d W / d x
    ! # dwdy = d W / d y
    ! # dwdz = d W / d z
    ! ###################################################################
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), u, dudx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), u, dudy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), u, dudz)

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), v, dvdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), v, dvdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), v, dvdz)

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, g(1), w, dwdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, g(2), w, dwdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, g(3), w, dwdz)

    ! ###################################################################
    ! Vorticity
    ! ###################################################################
    p_wrk3d = dwdy - dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortx(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = (p_wrk3d(:, :, k) - vortx(k))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortx2(1), wrk1d)

    p_wrk3d = dudz - dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vorty(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = (p_wrk3d(:, :, k) - vorty(k))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vorty2(1), wrk1d)

    p_wrk3d = dvdx - dudy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortz(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = (p_wrk3d(:, :, k) - vortz(k))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortz2(1), wrk1d)

    ! ###################################################################
    ! Derivatives Fluctuations
    ! ###################################################################
    ! -------------------------------------------------------------------
    ! Longitudinal terms
    p_wrk3d = dudx*dudx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_x2(1), wrk1d)
    p_wrk3d = p_wrk3d*dudx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_x3(1), wrk1d)
    p_wrk3d = p_wrk3d*dudx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_x4(1), wrk1d)

    p_wrk3d = dvdy*dvdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_y2(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_y3(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_y4(1), wrk1d)

    do k = 1, kmax
        p_wrk3d(:, :, k) = (dwdz(:, :, k) - rW_z(k))*(dwdz(:, :, k) - rW_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_z2(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = p_wrk3d(:, :, k)*(dwdz(:, :, k) - rW_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_z3(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = p_wrk3d(:, :, k)*(dwdz(:, :, k) - rW_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_z4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Lateral terms U
    p_wrk3d = dudy*dudy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_y2(1), wrk1d)
    p_wrk3d = p_wrk3d*dudy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_y3(1), wrk1d)
    p_wrk3d = p_wrk3d*dudy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_y4(1), wrk1d)

    do k = 1, kmax
        p_wrk3d(:, :, k) = (dudz(:, :, k) - rU_z(k))*(dudz(:, :, k) - rU_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_z2(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = p_wrk3d(:, :, k)*(dudz(:, :, k) - rU_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_z3(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = p_wrk3d(:, :, k)*(dudz(:, :, k) - rU_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_z4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Lateral terms V
    p_wrk3d = dvdx*dvdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_x2(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_x3(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_x4(1), wrk1d)

    do k = 1, kmax
        p_wrk3d(:, :, k) = (dvdz(:, :, k) - rU_z(k))*(dvdz(:, :, k) - rV_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_z2(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = p_wrk3d(:, :, k)*(dvdz(:, :, k) - rV_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_z3(1), wrk1d)
    do k = 1, kmax
        p_wrk3d(:, :, k) = p_wrk3d(:, :, k)*(dvdz(:, :, k) - rV_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_z4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Lateral terms W
    p_wrk3d = dwdx*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_x2(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_x3(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_x4(1), wrk1d)

    p_wrk3d = dwdy*dwdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_y2(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_y3(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_y4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Dilatation fluctuation
    p_wrk3d = dudx + dvdy + dwdz
    do k = 1, kmax
        p_wrk3d(:, :, k) = (p_wrk3d(:, :, k) - rW_z(k))*(p_wrk3d(:, :, k) - rW_z(k))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_ii2(1), wrk1d)

    ! ##################################################################
    ! Mean viscous dissipation rate
    ! ##################################################################
    p_wrk3d = dudx**2 + dvdy**2 + dwdz**2 - (dudx + dvdy + dwdz)**2/3.0_wp + &
              0.5_wp*((dudy + dvdx)**2 + (dudz + dwdx)**2 + (dvdz + dwdy)**2)

    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Phi(1), wrk1d)
    Phi(:) = Phi(:)*visc*2.0_wp

    ! ###################################################################
    ! Dissipation Terms; final operation after viscous terms below
    ! ###################################################################
    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dudx + (dudy + dvdx)*dudy + (dudz + dwdx)*dudz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Exx(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dvdy*2.0_wp - p_wrk3d)*dvdy + (dudy + dvdx)*dvdx + (dvdz + dwdy)*dvdz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Eyy(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dwdz*2.0_wp - p_wrk3d)*dwdz + (dwdy + dvdz)*dwdy + (dwdx + dudz)*dwdx
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Ezz(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dvdx + (dudy + dvdx)*dvdy + (dudz + dwdx)*dvdz &
              + (dvdy*2.0_wp - p_wrk3d)*dudy + (dudy + dvdx)*dudx + (dvdz + dwdy)*dudz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Exy(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dwdx + (dudy + dvdx)*dwdy + (dudz + dwdx)*dwdz &
              + (dwdz*2.0_wp - p_wrk3d)*dudz + (dudz + dwdx)*dudx + (dvdz + dwdy)*dudy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Exz(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dvdy*2.0_wp - p_wrk3d)*dwdy + (dudy + dvdx)*dwdx + (dvdz + dwdy)*dwdz &
              + (dwdz*2.0_wp - p_wrk3d)*dvdz + (dudz + dwdx)*dvdx + (dvdz + dwdy)*dvdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Eyz(1), wrk1d)

    ! ##################################################################
    ! Viscous shear-stress tensor
    ! ##################################################################
    p_wrk3d = dudz + dwdx
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_xz(1), wrk1d)
    do k = 1, kmax
        dudz(:, :, k) = p_wrk3d(:, :, k) - Tau_xz(k)                ! fluctuation tau13'
    end do
    Tau_xz(:) = Tau_xz(:)*visc

    p_wrk3d = dvdz + dwdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yz(1), wrk1d)
    do k = 1, kmax
        dvdz(:, :, k) = p_wrk3d(:, :, k) - Tau_yz(k)                ! fluctuation tau23'
    end do
    Tau_yz(:) = Tau_yz(:)*visc

    p_wrk3d = dwdz*2.0_wp - dudx - dvdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_zz(1), wrk1d)
    do k = 1, kmax
        dwdz(:, :, k) = (p_wrk3d(:, :, k) - Tau_zz(k))*c23          ! fluctuation tau33'
    end do
    Tau_zz(:) = Tau_zz(:)*visc*c23

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tau_xz(1), Tau_xz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tau_yz(1), Tau_yz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tau_zz(1), Tau_zz_z(1))

    ! -------------------------------------------------------------------
    ! Contribution to turbulent transport terms
    do k = 1, kmax
        p_wrk3d(:, :, k) = dudz(:, :, k)*(u(:, :, k) - fU(k))   ! -2*u'*tau13'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(:), wrk1d)
    Txxz(:) = Txxz(:) - aux(:)*visc*2.0_wp
    Tz3(:) = Tz3(:) - aux(:)*visc

    do k = 1, kmax
        p_wrk3d(:, :, k) = dvdz(:, :, k)*(v(:, :, k) - fV(k))   ! -2*v'*tau23'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(:), wrk1d)
    Tyyz(:) = Tyyz(:) - aux(:)*visc*2.0_wp
    Tz3(:) = Tz3(:) - aux(:)*visc

    do k = 1, kmax
        p_wrk3d(:, :, k) = dwdz(:, :, k)*(w(:, :, k) - fW(k))   ! -2*w'*tau33'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(:), wrk1d)
    Tzzz(:) = Tzzz(:) - aux(:)*visc*2.0_wp
    Tz3(:) = Tz3(:) - aux(:)*visc

    do k = 1, kmax
        p_wrk3d(:, :, k) = dvdz(:, :, k)*(u(:, :, k) - fU(k)) + dudz(:, :, k)*(v(:, :, k) - fV(k))! -u'*tau23' -v'*tau13'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(:), wrk1d)
    Txyz(:) = Txyz(:) - aux(:)*visc

    do k = 1, kmax
        p_wrk3d(:, :, k) = dwdz(:, :, k)*(u(:, :, k) - fU(k)) + dudz(:, :, k)*(w(:, :, k) - fW(k))! -u'*tau33' -w'*tau13'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(:), wrk1d)
    Txzz(:) = Txzz(:) - aux(:)*visc

    do k = 1, kmax
        p_wrk3d(:, :, k) = dwdz(:, :, k)*(v(:, :, k) - fV(k)) + dwdz(:, :, k)*(w(:, :, k) - fW(k))! -v'*tau23' -w'*tau33'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(:), wrk1d)
    Tyzz(:) = Tyzz(:) - aux(:)*visc

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Txxz(1), Txxz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tyyz(1), Tyyz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tzzz(1), Tzzz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Txyz(1), Txyz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Txzz(1), Txzz_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tyzz(1), Tyzz_z(1))

    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tz1(1), Tz1_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tz2(1), Tz2_z(1))
    call OPR_Partial_Z(OPR_P1, 1, 1, kmax, g(3), Tz3(1), Tz3_z(1))

    ! -------------------------------------------------------------------
    ! Contribution to dissipation
    Exx(:) = (Exx(:)*visc - Tau_xz(:)*rU_z(:))*2.0_wp
    Eyy(:) = (Eyy(:)*visc - Tau_yz(:)*rV_z(:))*2.0_wp
    Ezz(:) = (Ezz(:)*visc - Tau_zz(:)*rW_z(:))*2.0_wp
    Exy(:) = Exy(:)*visc - Tau_xz(:)*rV_z(:) - Tau_yz(:)*rU_z(:)
    Exz(:) = Exz(:)*visc - Tau_xz(:)*rW_z(:) - Tau_zz(:)*rU_z(:)
    Eyz(:) = Eyz(:)*visc - Tau_yz(:)*rW_z(:) - Tau_zz(:)*rV_z(:)

    ! ###################################################################
    ! Complete budget equations
    ! ###################################################################
    ! Rij Convective Terms
    Cxx(:) = -fV(:)*Rxx_z(:)
    Cyy(:) = -fV(:)*Ryy_z(:)
    Czz(:) = -fV(:)*Rzz_z(:)
    Cxy(:) = -fV(:)*Rxy_z(:)
    Cxz(:) = -fV(:)*Rxz_z(:)
    Cyz(:) = -fV(:)*Ryz_z(:)

    ! Rij Production Terms
    Pxx(:) = -2.0_wp*Rxz(:)*fU_z(:)
    Pyy(:) = -2.0_wp*Ryz(:)*fV_z(:)
    Pzz(:) = -2.0_wp*Rzz(:)*fW_z(:)
    Pxy(:) = -(Rxz(:)*fV_z(:) + Ryz(:)*fU_z(:))
    Pxz(:) = -(Rxz(:)*fW_z(:) + Rzz(:)*fU_z(:))
    Pyz(:) = -(Ryz(:)*fW_z(:) + Rzz(:)*fV_z(:))

    ! Rij Pressure Variable-Density Terms
    Gxx(:) = 0.0_wp
    Gyy(:) = 0.0_wp
    Gzz(:) = 2.0_wp*rWf(:)*rP_z(:)
    Gxy(:) = 0.0_wp
    Gxz(:) = rUf(:)*rP_z(:)
    Gyz(:) = rVf(:)*rP_z(:)

    ! Rij Viscous Variable-Density Terms
    Dxx(:) = 2.0_wp*rUf(:)*Tau_xz_z(:)
    Dyy(:) = 2.0_wp*rVf(:)*Tau_yz_z(:)
    Dzz(:) = 2.0_wp*rWf(:)*Tau_zz_z(:)
    Dxy(:) = rUf(:)*Tau_yz_z(:) + rVf(:)*Tau_xz_z(:)
    Dxz(:) = rUf(:)*Tau_zz_z(:) + rWf(:)*Tau_xz_z(:)
    Dyz(:) = rVf(:)*Tau_zz_z(:) + rWf(:)*Tau_yz_z(:)

    ! ! Rij Coriolis Terms
    ! if (coriolis%active(1) .and. coriolis%active(3)) then ! contribution from angular velocity Oy
    !     dummy = coriolis%vector(2)
    !     Fxx(:) = dummy*2.0_wp*Rxz(:)
    !     Fyy(:) = 0.0_wp
    !     Fzz(:) = -dummy*2.0_wp*Rxz(:)
    !     Fxy(:) = dummy*Ryz(:)
    !     Fxz(:) = dummy*(Rzz(:) - Rxx(:))
    !     Fyz(:) = -dummy*Rxy(:)
    ! end if

    ! Rij Transient terms
    select case (nse_eqns)

    case (DNS_EQNS_BOUSSINESQ)
        aux(:) = 1.0_wp

    case (DNS_EQNS_ANELASTIC)
        aux(:) = 1.0_wp/rbackground(:)

    case (DNS_EQNS_COMPRESSIBLE)
        aux(:) = 1.0_wp/rR(:)

    end select

    Rxx_t(:) = -Fxx(:) + Bxx(:) + Cxx(:) + Pxx(:) - Exx(:) + (PIxx(:) - Txxz_z(:) - Gxx(:) + Dxx(:))*aux(:)
    Ryy_t(:) = -Fyy(:) + Byy(:) + Cyy(:) + Pyy(:) - Eyy(:) + (PIyy(:) - Tyyz_z(:) - Gyy(:) + Dyy(:))*aux(:)
    Rzz_t(:) = -Fzz(:) + Bzz(:) + Czz(:) + Pzz(:) - Ezz(:) + (PIzz(:) - Tzzz_z(:) - Gzz(:) + Dzz(:))*aux(:)
    Rxy_t(:) = -Fxy(:) + Bxy(:) + Cxy(:) + Pxy(:) - Exy(:) + (PIxy(:) - Txyz_z(:) - Gxy(:) + Dxy(:))*aux(:)
    Rxz_t(:) = -Fxz(:) + Bxz(:) + Cxz(:) + Pxz(:) - Exz(:) + (PIxz(:) - Txzz_z(:) - Gxz(:) + Dxz(:))*aux(:)
    Ryz_t(:) = -Fyz(:) + Byz(:) + Cyz(:) + Pyz(:) - Eyz(:) + (PIyz(:) - Tyzz_z(:) - Gyz(:) + Dyz(:))*aux(:)

    ! Kinetic energy equation
    Tke(:) = 0.5_wp*(Rxx(:) + Ryy(:) + Rzz(:))

    Buo(:) = 0.5_wp*(Bxx(:) + Byy(:) + Bzz(:))
    Con(:) = 0.5_wp*(Cxx(:) + Cyy(:) + Czz(:))
    Prd(:) = 0.5_wp*(Pxx(:) + Pyy(:) + Pzz(:))
    Pi(:) = 0.5_wp*(PIxx(:) + PIyy(:) + PIzz(:))
    Eps(:) = 0.5_wp*(Exx(:) + Eyy(:) + Ezz(:))
    Tz_z(:) = 0.5_wp*(Txxz_z(:) + Tyyz_z(:) + Tzzz_z(:))
    Gkin(:) = 0.5_wp*(Gxx(:) + Gyy(:) + Gzz(:))
    Dkin(:) = 0.5_wp*(Dxx(:) + Dyy(:) + Dzz(:))

    Tke_t(:) = Buo(:) + Con(:) + Prd(:) - Eps(:) + (-Tz_z(:) + Pi(:) - Gkin(:) + Dkin(:))*aux(:)

    ! ###################################################################
    ! Output
    ! ###################################################################
    ! 14 t-dependent variables, for consistency with old format
    ! ng = ng +1
    ! groupname(ng) = ''
    ! varname(ng)   = 'dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy'
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''

    write (name, *) itime; name = 'avg'//trim(adjustl(name))
    call IO_WRITE_AVERAGES(name, itime, rtime, kmax, nv, ng, z%nodes, varname, groupname, mean2d)

    return
end subroutine AVG_FLOW_XZ
