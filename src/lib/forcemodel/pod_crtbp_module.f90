module pod_crtbp_module
    use pod_global, only: DP
    use pod_dace_classes
    implicit none
    private

    real(DP), public :: crtbp_mu = 0.01_DP
    public :: set_crtbp_mu
    public :: crtbp_derivatives_real
    public :: crtbp_derivatives_da
    public :: CrtbpForceTempPool

    type, public :: CrtbpForceTempPool
        type(DA) :: r1, r2, v1, v2
        type(DA) :: Omega_x, Omega_y, Omega_z
        type(DA) :: tmp1, tmp2, tmp3
    contains
        procedure :: init   => crtbp_force_pool_init
        procedure :: destroy => crtbp_force_pool_destroy
    end type CrtbpForceTempPool

contains

    subroutine set_crtbp_mu(mu)
        real(DP), intent(in) :: mu
        crtbp_mu = mu
    end subroutine set_crtbp_mu

    subroutine crtbp_force_pool_init(this)
        class(CrtbpForceTempPool), intent(inout) :: this
        call this%r1%init();      call this%r2%init()
        call this%v1%init();      call this%v2%init()
        call this%Omega_x%init(); call this%Omega_y%init(); call this%Omega_z%init()
        call this%tmp1%init();    call this%tmp2%init();    call this%tmp3%init()
    end subroutine crtbp_force_pool_init

    subroutine crtbp_force_pool_destroy(this)
        class(CrtbpForceTempPool), intent(inout) :: this
        call this%r1%destroy();      call this%r2%destroy()
        call this%v1%destroy();      call this%v2%destroy()
        call this%Omega_x%destroy(); call this%Omega_y%destroy(); call this%Omega_z%destroy()
        call this%tmp1%destroy();    call this%tmp2%destroy();    call this%tmp3%destroy()
    end subroutine crtbp_force_pool_destroy

    subroutine crtbp_derivatives_real(x, dxdt, t)
        real(DP), dimension(6), intent(in)  :: x
        real(DP), dimension(6), intent(out) :: dxdt
        real(DP), intent(in) :: t

        real(DP) :: r1, r2, v1, v2
        real(DP) :: Omega_x, Omega_y, Omega_z
        real(DP) :: mu

        mu = crtbp_mu

        r1 = sqrt((x(1) + mu)**2 + x(2)**2 + x(3)**2)
        r2 = sqrt((x(1) - 1.0_DP + mu)**2 + x(2)**2 + x(3)**2)
        v1 = (1.0_DP - mu) / (r1**3)
        v2 = mu / (r2**3)

        Omega_x = x(1) - (x(1) + mu) * v1 - v2 * (x(1) - 1.0_DP + mu)
        Omega_y = x(2) * (1.0_DP - v1 - v2)
        Omega_z = -x(3) * (v1 + v2)

        dxdt(1) = x(4)
        dxdt(2) = x(5)
        dxdt(3) = x(6)
        dxdt(4) = Omega_x + 2.0_DP * x(5)
        dxdt(5) = Omega_y - 2.0_DP * x(4)
        dxdt(6) = Omega_z
    end subroutine crtbp_derivatives_real

    subroutine crtbp_derivatives_da(x, dxdt, t, pool)
        type(DA), dimension(6), intent(in)    :: x
        type(DA), dimension(6), intent(inout) :: dxdt
        real(DP), intent(in) :: t
        type(CrtbpForceTempPool), intent(inout) :: pool

        real(DP) :: mu

        mu = crtbp_mu

        ! r1 = sqrt((x1+mu)^2 + x2^2 + x3^2)
        pool%tmp1 = x(1) + mu
        pool%r1 = pool%tmp1 * pool%tmp1 + x(2) * x(2) + x(3) * x(3)
        call da_sqrt_sub(pool%r1, pool%r1)

        ! r2 = sqrt((x1-1+mu)^2 + x2^2 + x3^2)
        pool%tmp2 = x(1) - 1.0_DP + mu
        pool%r2 = pool%tmp2 * pool%tmp2 + x(2) * x(2) + x(3) * x(3)
        call da_sqrt_sub(pool%r2, pool%r2)

        ! v1 = (1-mu) / r1^3
        pool%tmp1 = pool%r1 * pool%r1 * pool%r1
        pool%v1 = (1.0_DP - mu) / pool%tmp1

        ! v2 = mu / r2^3
        pool%tmp1 = pool%r2 * pool%r2 * pool%r2
        pool%v2 = mu / pool%tmp1

        ! Omega_x = x1 - (x1+mu)*v1 - v2*(x1-1+mu)
        pool%Omega_x = x(1) - (x(1) + mu) * pool%v1 - pool%v2 * (x(1) - 1.0_DP + mu)

        ! Omega_y = x2 * (1 - v1 - v2)
        pool%tmp1 = 1.0_DP - pool%v1 - pool%v2
        pool%Omega_y = x(2) * pool%tmp1

        ! Omega_z = -x3 * (v1 + v2)
        pool%tmp1 = pool%v1 + pool%v2
        pool%Omega_z = -x(3) * pool%tmp1

        dxdt(1) = x(4)
        dxdt(2) = x(5)
        dxdt(3) = x(6)
        dxdt(4) = pool%Omega_x + 2.0_DP * x(5)
        dxdt(5) = pool%Omega_y - 2.0_DP * x(4)
        dxdt(6) = pool%Omega_z
    end subroutine crtbp_derivatives_da

end module pod_crtbp_module
