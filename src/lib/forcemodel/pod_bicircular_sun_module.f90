! Deterministic Sun geometry in the Earth-Moon barycentric rotating frame.
! This module supplies NO gravitational acceleration.
module pod_bicircular_sun_module
    use pod_global, only: DP
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: BicircularSunConfig, bicircular_sun_position, bicircular_sun_phase
    public :: sun_elapsed_tau, validate_bicircular_sun

    type :: BicircularSunConfig
        ! Match app/run_CRTBP_uprop.f90; its LU is km, this API is strictly SI.
        real(DP) :: du_m = 384400000.0_DP
        real(DP) :: tu_s = 375190.26189464843_DP
        real(DP) :: au_m = 149597870700.0_DP
        real(DP) :: radius_m = 149597870700.0_DP
        real(DP) :: period_s = 31558149.8_DP
        real(DP) :: alpha0_rad = 0.0_DP ! phase at tau=0, i.e. physical epoch t0
    end type
contains
    subroutine validate_bicircular_sun(cfg)
        type(BicircularSunConfig), intent(in) :: cfg
        real(DP) :: positive(5)
        positive = [cfg%du_m,cfg%tu_s,cfg%au_m,cfg%radius_m,cfg%period_s]
        if (.not.all(ieee_is_finite(positive)) .or. any(positive<=0.0_DP)) &
            error stop 'bicircular Sun: DU, TU, AU, radius and period must be finite and positive'
        if (.not.ieee_is_finite(cfg%alpha0_rad)) error stop 'bicircular Sun: non-finite initial phase'
    end subroutine

    function sun_elapsed_tau(t_seconds,t0_seconds,cfg) result(tau)
        real(DP), intent(in) :: t_seconds,t0_seconds
        type(BicircularSunConfig), intent(in) :: cfg
        real(DP) :: tau
        call validate_bicircular_sun(cfg)
        if (.not.all(ieee_is_finite([t_seconds,t0_seconds]))) error stop 'bicircular Sun: non-finite epoch'
        tau = (t_seconds-t0_seconds)/cfg%tu_s
    end function

    function bicircular_sun_phase(tau,cfg) result(alpha)
        real(DP), intent(in) :: tau
        type(BicircularSunConfig), intent(in) :: cfg
        real(DP) :: alpha
        call validate_bicircular_sun(cfg)
        if (.not.ieee_is_finite(tau)) error stop 'bicircular Sun: non-finite tau'
        ! n=1/TU. The inertial Sun angle minus the rotating-frame angle.
        alpha = cfg%alpha0_rad + (2.0_DP*acos(-1.0_DP)*cfg%tu_s/cfg%period_s-1.0_DP)*tau
    end function

    subroutine bicircular_sun_position(tau,cfg,r_sun_nd)
        real(DP), intent(in) :: tau
        type(BicircularSunConfig), intent(in) :: cfg
        real(DP), intent(out) :: r_sun_nd(3)
        real(DP) :: alpha
        alpha = bicircular_sun_phase(tau,cfg)
        r_sun_nd = (cfg%radius_m/cfg%du_m)*[cos(alpha),sin(alpha),0.0_DP]
    end subroutine
end module pod_bicircular_sun_module
