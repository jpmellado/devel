subroutine Tlab_Diagnostic(nx, ny, nz, s)
    use TLab_Constants, only: wp, wi
    use Thermodynamics, only: imode_thermo, THERMO_TYPE_ANELASTIC, THERMO_TYPE_COMPRESSIBLE
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_AirWater, only: inb_scal_T
    use Thermo_Anelastic
    use Microphysics
    implicit none

    integer(wi), intent(in) :: nx, ny, nz
    real(wp), intent(inout) :: s(nx*ny*nz, *)

    !########################################################################
    select case (imode_thermo)
    case (THERMO_TYPE_ANELASTIC)
        if (imixture == MIXT_TYPE_AIRWATER .and. evaporationProps%type == TYPE_EVA_EQUILIBRIUM) then
            call Thermo_Anelastic_EquilibriumPH(nx, ny, nz, s(:, 2), s(:, 1))
        end if
        call Thermo_Anelastic_T(nx, ny, nz, s, s(:, inb_scal_T))

    case (THERMO_TYPE_COMPRESSIBLE)

    end select

    return
end subroutine Tlab_Diagnostic
