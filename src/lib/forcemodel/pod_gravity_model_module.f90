!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!> # POD Gravity Model Module
!> 
!> This module provides high-order spherical harmonics gravity field modeling
!> for the POD Fortran precision orbit determination (POD) system. It calculates
!> the gravitational acceleration induced by the non-spherical mass distribution
!> of a central body (e.g., Earth, Moon).
!> 
!> ## Features
!> 
!> - **Model Parsing**: Reads standard gravity field files (e.g., GGM05C.GEO).
!> - **Zonal Harmonics**: Computes acceleration from zonal terms (m=0) via `f_zonal`.
!> - **Tesseral Harmonics**: Computes acceleration from sectorial and tesseral terms (m>0) via `f_tesseral`.
!> - **Legendre Polynomials**: Utilizes fully normalized associated Legendre polynomials
!>   with highly stable recurrence relations.
!> 
!> ## Dependencies
!> 
!> - `pod_global`: For precision definitions (`DP`).
!> - `pod_dace_classes`: For DA vector types used in high-order uncertainty propagation.
!> - `pod_spice`: For time conversions and ephemeris data access (if needed for gravity field parameters).
!> 
!> ## Input Files
!> 
!> - Gravity field model files loPODed in the `data/` directory:
!>   - `GGM05C.GEO` (Earth gravity model)
!>   - `gggrx_0660pm_sha.tab` (Moon gravity model)
!> 
!> ## Authors
!> 
!> - **Original Author**: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> - **Refactored by**: Wenxuan (Integrated into POD POD architecture)
!> 
!> ## References
!> 
!> - Liu Lin, "Satellite Orbit Mechanics and AppliPODions" (《卫星轨道力学与应用》).
!> - Base algorithm logic referenced from Hu Shoucun.
!> 
!> @note This module currently evaluates accelerations using pure double precision (`real(DP)`).
!>       For high-order uncertainty mapping, a Differential Algebra (DA) overloaded version
!>       should be implemented alongside this baseline.
!> 
!> @warning The `is_gravity_model_loaded` flag in the force model wrapper must be 
!>          checked before calling evaluation routines to prevent uninitialized memory access.
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------

module pod_gravity_model_module
    use pod_global, only: DP
    use pod_dace_classes
    implicit none
    
    real(DP) :: g_m, r_m, g_e, r_e, gmm, rm
    integer, parameter :: ndeg_max = 100
    
    real(DP), dimension(ndeg_max) :: cl0
    real(DP), dimension(100), public :: cl0e
    real(DP), dimension(100), public :: cl0m
    real(DP), dimension(ndeg_max, ndeg_max) :: clm, slm
    real(DP), dimension(ndeg_max, ndeg_max) :: clme, clmm, slme, slmm
    
    type, public :: gravity_field
        real(DP), dimension(3) :: dr
        type(AlgebraicVector) :: dr_da
        integer :: ncs, cen_body
    contains
        procedure, public, pass(gf) :: read_gravity_field
        procedure, public, pass(gf) :: f_zonal, f_zonal_da
        procedure, public, pass(gf) :: f_tesseral, f_tesseral_da
    end type gravity_field
    
    ! =========================================================
    ! 临时变量池：集中管理 DA 运算所需的临时变量，
    ! 避免在循环中反复申请/释放 DA 句柄导致内存泄露。
    ! =========================================================
    type :: GravityTempPool
        ! --- DA 标量临时变量 ---
        type(DA) :: r, r3, rl3, zr, rg2, ratio_tmp
        type(DA) :: tmp2, tmp3
        type(DA) :: eta, coslg, sinlg, w11_da, w21_da, dplm
        type(DA) :: sqrt_res
        type(DA) :: tmp_da1, tmp_da2, tmp_da3, tmp_da4, tmp_da5
        type(DA) :: tmp_da6, tmp_da7, tmp_da8
        ! --- AlgebraicVector 临时变量 ---
        type(AlgebraicVector) :: h, k_v, g_v, w1, w2, dr_da_tmp
        type(AlgebraicVector) :: tmp_vec1, tmp_vec2, tmp_vec3, tmp_vec4
        type(AlgebraicVector) :: tmp_vec5, tmp_vec6, tmp_vec7
    contains
        procedure :: init => pool_init
        procedure :: destroy => pool_destroy
    end type GravityTempPool

contains

    ! =========================================================
    ! 临时变量池初始化
    ! =========================================================
    subroutine pool_init(this, n)
        class(GravityTempPool), intent(inout) :: this
        integer, intent(in) :: n
        call this%r%init(); call this%r3%init(); call this%rl3%init()
        call this%zr%init(); call this%rg2%init(); call this%ratio_tmp%init()
        call this%tmp2%init(); call this%tmp3%init()
        call this%eta%init(); call this%coslg%init(); call this%sinlg%init()
        call this%w11_da%init(); call this%w21_da%init(); call this%dplm%init()
        call this%sqrt_res%init()
        call this%tmp_da1%init();  call this%tmp_da2%init();  call this%tmp_da3%init()
        call this%tmp_da4%init();  call this%tmp_da5%init();  call this%tmp_da6%init()
        call this%tmp_da7%init();  call this%tmp_da8%init()
        call this%h%init(n); call this%k_v%init(n); call this%g_v%init(n)
        call this%w1%init(n); call this%w2%init(n); call this%dr_da_tmp%init(n)
        call this%tmp_vec1%init(n); call this%tmp_vec2%init(n)
        call this%tmp_vec3%init(n); call this%tmp_vec4%init(n)
        call this%tmp_vec5%init(n); call this%tmp_vec6%init(n)
        call this%tmp_vec7%init(n)
    end subroutine pool_init

    ! =========================================================
    ! 临时变量池销毁
    ! =========================================================
    subroutine pool_destroy(this)
        class(GravityTempPool), intent(inout) :: this
        call this%r%destroy(); call this%r3%destroy(); call this%rl3%destroy()
        call this%zr%destroy(); call this%rg2%destroy(); call this%ratio_tmp%destroy()
        call this%tmp2%destroy(); call this%tmp3%destroy()
        call this%eta%destroy(); call this%coslg%destroy(); call this%sinlg%destroy()
        call this%w11_da%destroy(); call this%w21_da%destroy(); call this%dplm%destroy()
        call this%sqrt_res%destroy()
        call this%tmp_da1%destroy();  call this%tmp_da2%destroy();  call this%tmp_da3%destroy()
        call this%tmp_da4%destroy();  call this%tmp_da5%destroy();  call this%tmp_da6%destroy()
        call this%tmp_da7%destroy();  call this%tmp_da8%destroy()
        call this%h%destroy(); call this%k_v%destroy(); call this%g_v%destroy()
        call this%w1%destroy(); call this%w2%destroy(); call this%dr_da_tmp%destroy()
        call this%tmp_vec1%destroy(); call this%tmp_vec2%destroy()
        call this%tmp_vec3%destroy(); call this%tmp_vec4%destroy()
        call this%tmp_vec5%destroy(); call this%tmp_vec6%destroy()
        call this%tmp_vec7%destroy()
    end subroutine pool_destroy

    ! =========================================================
    ! 读取引力场模型文件
    ! =========================================================
    subroutine read_gravity_field(gf)
        class(gravity_field), intent(in) :: gf
        character(100) :: filename
        integer :: i1, i2, nstat, ios
        real(DP), dimension(2) :: tmp
        
        if (gf%cen_body == 3)  filename = './kernels/gravity_models/GGM05C.GEO'
        if (gf%cen_body == 10) filename = './kernels/gravity_models/gggrx_0660pm_sha.tab'  
        
        cl0 = 0.0_DP; clm = 0.0_DP; slm = 0.0_DP
        
        open(101, file=trim(filename), iostat=ios)
        if (ios == 0) then
            if (gf%cen_body == 3) then
                read(101, *) g_e, r_e
                g_e = g_e / 1e9_DP
                r_e = r_e / 1e3_DP
            end if
            if (gf%cen_body == 10) then
                read(101, *) g_m, r_m
                g_m = g_m / 1e9_DP
                r_m = r_m / 1e3_DP
            end if
        else
            write(*,*) '错误: 无法打开引力场模型文件: ', trim(filename)
            stop
        end if
        
        do while (.true.)
            read(101, *, iostat=nstat) i1, i2, tmp
            if (nstat /= 0) exit
            if (i1 > gf%ncs .or. i2 > gf%ncs) exit
            if (i1 > ndeg_max .or. i2 > ndeg_max) exit
            
            if (i2 == 0) then
                if (gf%cen_body == 3) cl0e(i1) = tmp(1)
                if (gf%cen_body == 10) cl0m(i1) = tmp(1)
            else 
                if (gf%cen_body == 3) then
                    clme(i1,i2) = tmp(1); slme(i1,i2) = tmp(2)
                else
                    clmm(i1,i2) = tmp(1); slmm(i1,i2) = tmp(2)
                end if
            end if
        end do    
        close(101)
    end subroutine read_gravity_field

    ! =========================================================
    ! 实数版本：带谐项加速度计算
    ! =========================================================
    subroutine f_zonal(gf, fl)
        class(gravity_field), intent(in) :: gf
        integer :: l
        real(DP), dimension(gf%ncs) :: pl
        real(DP) :: r, r3, rl3, zr, rg2, tmp1, tmp2, tmp3
        real(DP), dimension(3) :: h, fl
        
        if (gf%cen_body == 3) then
            cl0 = cl0e; rm = r_e; gmm = g_e
        end if
        if (gf%cen_body == 10) then
            cl0 = cl0m; rm = r_m; gmm = g_m
        end if
        
        fl = 0.0_DP
        r = norm2(gf%dr)
        zr = gf%dr(3) / r
        r3 = r**3_DP
        rg2 = gf%dr(1)**2_DP + gf%dr(2)**2_DP
        h = (/gf%dr(1)*gf%dr(3), gf%dr(2)*gf%dr(3), -rg2/)

        if (gf%ncs < 2) return
        call plx(gf%ncs, zr, pl)
        
        do l = 2, gf%ncs
            rl3 = (r/rm)**real(l, DP) * r3
            tmp1 = sqrt((2.0_DP*l+1.0_DP)/(2.0_DP*real(l, DP)-1.0_DP))
            tmp2 = (l+1.0_DP) * pl(l)
            
            if (rg2/r > 1.0e-14_DP) then 
                tmp3 = real(l, DP)*r/rg2 * (tmp1*pl(l-1) - zr*pl(l))
            else 
                tmp3 = 1.0_DP
            end if
            fl = fl - cl0(l) * (tmp2*gf%dr + tmp3*h) / rl3
        end do
        fl = fl * gmm 
    end subroutine f_zonal
   
    ! =========================================================
    ! 实数版本：田谐项加速度计算
    ! =========================================================
    subroutine f_tesseral(gf, flm)
        class(gravity_field), intent(in) :: gf
        integer :: l, m
        real(DP), dimension(gf%ncs) :: cosmlg, sinmlg
        real(DP), dimension(gf%ncs, gf%ncs+1) :: plm
        real(DP) :: coslg, sinlg, r, r2, zr, w11, w21, eta, dplm
        real(DP), dimension(3) :: k_v, g_v, w1, w2, flm
        
        if (gf%ncs < 1) return
        
        if (gf%cen_body == 3) then
            clm = clme; slm = slme
        else 
            clm = clmm; slm = slmm
        end if
            
        flm = 0.0_DP
        r = norm2(gf%dr)
        r2 = sqrt(gf%dr(1)**2_DP + gf%dr(2)**2_DP)
        zr = gf%dr(3) / r
        coslg = gf%dr(1) / r2
        sinlg = gf%dr(2) / r2
        cosmlg(1) = coslg; sinmlg(1) = sinlg
    
        if (gf%ncs < 2) return
        if (gf%ncs > ndeg_max) stop
        
        call plmx(gf%ncs, zr, plm)
        
        k_v = (/0.0_DP, 0.0_DP, 1.0_DP/)
        g_v = (/-sinlg, coslg, 0.0_DP/)
        eta = sqrt(1.0_DP - zr**2_DP)

        if (eta < 1e-7_DP) stop
    
        do l = 2, gf%ncs
            if (l == 2) then
                cosmlg(l) = 2.0_DP*coslg*cosmlg(l-1) - 1.0_DP
                sinmlg(l) = 2.0_DP*coslg*sinmlg(l-1)
            else
                cosmlg(l) = 2.0_DP*coslg*cosmlg(l-1) - cosmlg(l-2)
                sinmlg(l) = 2.0_DP*coslg*sinmlg(l-1) - sinmlg(l-2)
            end if
            
            do m = 1, l
                w11 = sqrt((l+m+1.0_DP)*(l-m))
                w21 = m*zr/eta
                dplm = (w11*plm(l,m+1) - w21*plm(l,m))/eta
                w1 = (rm/r)**l/r**3_DP * (((l+1.0_DP)*plm(l,m)+zr*dplm)*gf%dr - r*dplm*k_v)
                w11 = clm(l,m)*cosmlg(m) + slm(l,m)*sinmlg(m)
                w2 = m/r2/r*(rm/r)**l * plm(l,m)*g_v
                w21 = clm(l,m)*sinmlg(m) - slm(l,m)*cosmlg(m)
                flm = flm - (w1*w11 + w2*w21)
            end do
        end do
        flm = flm * gmm
    end subroutine f_tesseral
    
    ! =========================================================
    ! 实数版本：Legendre 多项式
    ! =========================================================
    subroutine plx(n, zr, pl)
        integer, intent(in) :: n
        real(DP), intent(in) :: zr
        real(DP), intent(out) :: pl(n)
        integer :: l
        real(DP) :: w1, w2, l2
        
        pl = 0.0_DP
        pl(1) = sqrt(3.0_DP) * zr
        pl(2) = sqrt(5.0_DP)/2.0_DP * (3.0_DP*zr*zr - 1.0_DP)
        if (n < 3) return
        do l = 3, n
            l2 = 2.0_DP * real(l, DP)
            w1 = sqrt((l2+1.0_DP)/(l2-1.0_DP))
            w2 = sqrt((l2-1.0_DP)/(l2-3.0_DP))
            pl(l) = w1 * ((2.0_DP-1.0_DP/real(l, DP))*zr*pl(l-1) - w2*(1.0_DP-1.0_DP/real(l, DP))*pl(l-2))
        end do
    end subroutine plx

    ! =========================================================
    ! 实数版本：连带 Legendre 多项式
    ! =========================================================
    subroutine PlmX(n, zr, plm)
        integer, intent(in) :: n
        real(DP), intent(in) :: zr
        real(DP), intent(out) :: plm(n,n)
        integer :: l, m
        real(DP) :: eta, w1, w2, l2
        
        eta = sqrt(1.0_DP - zr**2_DP)
        plm = 0.0_DP
        plm(1,1) = sqrt(3.0_DP) * eta
        plm(2,1) = sqrt(5.0_DP) * zr * plm(1,1)
        plm(2,2) = sqrt(5.0_DP)/2.0_DP * eta * plm(1,1)
        if (n < 3) return
        do l = 3, n
            l2 = 2.0_DP * real(l, DP)
            do m = 1, l-1
                w1 = (l2+1.0_DP)*(l2-1.0_DP) / (real(l+m, DP)*real(l-m, DP))
                w2 = (l2+1.0_DP)*(real(l-1+m, DP)*real(l-1-m, DP)) / ((l2-3.0_DP)*real(l+m, DP)*real(l-m, DP))
                plm(l,m) = sqrt(w1)*zr*plm(l-1,m) - sqrt(w2)*plm(l-2,m)
            end do
            plm(l,l) = sqrt((l2+1.0_DP)/l2) * eta * plm(l-1,l-1)
        end do
    end subroutine PlmX
    ! =========================================================
    ! DA 版本：带谐项加速度计算（显式运算，无临时变量泄露）
    ! =========================================================
    subroutine f_zonal_da(gf, fl)
        class(gravity_field), intent(in) :: gf
        type(AlgebraicVector), intent(inout) :: fl
        type(GravityTempPool) :: pool
        integer :: l
        type(DA), dimension(gf%ncs) :: pl
        real(DP) :: tmp1
        
        if (gf%cen_body == 3) then
            cl0 = cl0e; rm = r_e; gmm = g_e
        end if
        if (gf%cen_body == 10) then
            cl0 = cl0m; rm = r_m; gmm = g_m
        end if
        
        ! 初始化临时变量池
        call pool%init(3)
        
        call fl%init(3)
        call fl%sync_h()
        ! fl = 0.0_DP
        call vec_mul(0.0_DP, fl, fl)

        pool%dr_da_tmp = gf%dr_da !(复制)
        ! call vec_add(gf%dr_da, 0.0_DP*fl, pool%dr_da_tmp)
        ! call real_mul_vector_sub(1.0_DP, gf%dr_da, pool%dr_da_tmp)

        ! r = norm2(dr_da_tmp)  (使用 subroutine 版本避免临时句柄泄漏)
        call vector_norm2_sub(pool%dr_da_tmp, pool%r)

        ! zr = dr_da_tmp(3) / r
        call da_div(pool%dr_da_tmp%elements(3), pool%r, pool%zr)

        ! r3 = r**3  (使用 da_mul 链代替 da_pow_int_sub)
        call da_mul(pool%r, pool%r, pool%tmp_da1)
        call da_mul(pool%tmp_da1, pool%r, pool%r3)

        ! rg2 = dr_da_tmp(1)**2 + dr_da_tmp(2)**2
        call da_mul(pool%dr_da_tmp%elements(1), pool%dr_da_tmp%elements(1), pool%tmp_da1)
        call da_mul(pool%dr_da_tmp%elements(2), pool%dr_da_tmp%elements(2), pool%tmp_da2)
        call da_add(pool%tmp_da1, pool%tmp_da2, pool%rg2)

        ! h(1) = dr_da_tmp(1)*dr_da_tmp(3)
        call da_mul(pool%dr_da_tmp%elements(1), pool%dr_da_tmp%elements(3), pool%h%elements(1))
        ! h(2) = dr_da_tmp(2)*dr_da_tmp(3)
        call da_mul(pool%dr_da_tmp%elements(2), pool%dr_da_tmp%elements(3), pool%h%elements(2))
        ! h(3) = -rg2  (使用 unary_minus: 0 - rg2 不可用，改用 -1.0 * rg2)
        call da_mul(pool%rg2, -1.0_DP, pool%h%elements(3))

        if (gf%ncs < 2) then
            call pool%destroy()
            return
        end if

        ! 调用 plx_da 计算 pl 数组
        call plx_da(gf%ncs, pool%zr, pl, pool)
        
        do l = 2, gf%ncs
            ! rl3 = ((r/rm)**l) * r3
            ! 使用 da_pow_int_sub 避免临时句柄泄漏
            call da_div(pool%r, rm, pool%tmp_da1)
            call da_pow_int_sub(pool%tmp_da1, l, pool%tmp_da2)
            call da_mul(pool%tmp_da2, pool%r3, pool%rl3)

            tmp1 = sqrt((2.0_DP*l+1.0_DP)/(2.0_DP*real(l, DP)-1.0_DP))
            
            ! tmp2 = (l+1.0_DP) * pl(l)
            call da_mul(pl(l), (l+1.0_DP), pool%tmp2)
            
            ! ratio_tmp = rg2/r
            call da_div(pool%rg2, pool%r, pool%ratio_tmp)
            if (pool%ratio_tmp%cons() > 1.0e-14_DP) then 
                ! tmp3 = real(l, DP)*r/rg2 * (tmp1*pl(l-1) - zr*pl(l))
                call da_mul(pl(l-1), tmp1, pool%tmp_da3)
                call da_mul(pool%zr, pl(l), pool%tmp_da4)
                call da_sub(pool%tmp_da3, pool%tmp_da4, pool%tmp_da5)
                call da_div(pool%r, pool%rg2, pool%tmp_da6)
                call da_mul(pool%tmp_da6, real(l, DP), pool%tmp_da7)
                call da_mul(pool%tmp_da7, pool%tmp_da5, pool%tmp3)
            else 
                ! tmp3 = 1.0_DP
                call da_mul(pool%tmp3, 0.0_DP, pool%tmp3)
                call da_add(pool%tmp3, 1.0_DP, pool%tmp3)
            end if
            
            ! fl = fl - cl0(l) * (tmp2*dr_da_tmp + tmp3*h) / rl3
            ! 先计算 tmp2*dr_da_tmp
            call vec_mul(pool%tmp2, pool%dr_da_tmp, pool%tmp_vec1)
            ! 再计算 tmp3*h
            call vec_mul(pool%tmp3, pool%h, pool%tmp_vec2)
            ! 相加: tmp_vec1 + tmp_vec2
            call vec_add(pool%tmp_vec1, pool%tmp_vec2, pool%tmp_vec3)
            ! 除以 rl3
            call vec_div(pool%tmp_vec3, pool%rl3, pool%tmp_vec4)
            ! 乘以 cl0(l)
            call vec_mul(cl0(l), pool%tmp_vec4, pool%tmp_vec5)
            ! fl = fl - tmp_vec5
            call vec_sub(fl, pool%tmp_vec5, pool%tmp_vec6)
            ! call vec_add(pool%tmp_vec6, 0.0_DP*fl, fl)
            fl = pool%tmp_vec6
        end do
        
        ! fl = fl * gmm
        call vec_mul(gmm, fl, pool%tmp_vec1)
        ! call vec_add(pool%tmp_vec1, 0.0_DP*fl, fl)
        fl = pool%tmp_vec1
        
        ! 清理 pl 数组
        do l = 1, gf%ncs
            call pl(l)%destroy()
        end do
        
        ! 销毁临时变量池
        call pool%destroy()
    end subroutine f_zonal_da

    ! =========================================================
    ! DA 版本：田谐项加速度计算（显式运算，无临时变量泄露）
    ! =========================================================
    subroutine f_tesseral_da(gf, flm)
        class(gravity_field), intent(in) :: gf
        type(AlgebraicVector), intent(inout) :: flm
        type(GravityTempPool) :: pool
        integer :: l, m
        type(DA), dimension(gf%ncs) :: cosmlg, sinmlg
        type(DA), dimension(gf%ncs, gf%ncs+1) :: plm
        
        ! 在任何 return 之前，先建立安全防线！
        if (flm%size /= 3) call flm%init(3)
        
        ! 初始化临时变量池
        call pool%init(3)
        
        ! flm = 0.0_DP
        call flm%sync_h()
        call vec_mul(0.0_DP, flm, flm)
        
        ! --- 安全检查提前 ---
        if (gf%ncs < 1) then
            call pool%destroy()
            return
        end if
        
        if (gf%cen_body == 3) then
            clm = clme; slm = slme
        else 
            clm = clmm; slm = slmm
        end if

        ! dr_da_tmp = gf%dr_da
        ! call vec_add(gf%dr_da, 0.0_DP*flm, pool%dr_da_tmp)
        pool%dr_da_tmp = gf%dr_da
        
        ! r = norm2(dr_da_tmp)  (使用 subroutine 版本避免临时句柄泄漏)
        call vector_norm2_sub(pool%dr_da_tmp, pool%r)
        
        ! r2 = sqrt(dr_da_tmp(1)**2 + dr_da_tmp(2)**2)
        call da_mul(pool%dr_da_tmp%elements(1), pool%dr_da_tmp%elements(1), pool%tmp_da1)
        call da_mul(pool%dr_da_tmp%elements(2), pool%dr_da_tmp%elements(2), pool%tmp_da2)
        call da_add(pool%tmp_da1, pool%tmp_da2, pool%tmp_da3)
        ! sqrt 使用 subroutine 版本避免临时句柄泄漏
        call da_sqrt_sub(pool%tmp_da3, pool%sqrt_res)
        ! 用 pool%w11_da 暂存 r2
        call da_add(pool%sqrt_res, 0.0_DP, pool%w11_da)
        
        ! zr = dr_da_tmp(3) / r
        call da_div(pool%dr_da_tmp%elements(3), pool%r, pool%zr)
        
        ! coslg = dr_da_tmp(1) / r2
        call da_div(pool%dr_da_tmp%elements(1), pool%w11_da, pool%coslg)
        
        ! sinlg = dr_da_tmp(2) / r2
        call da_div(pool%dr_da_tmp%elements(2), pool%w11_da, pool%sinlg)

        ! cosmlg(1) = coslg; sinmlg(1) = sinlg
        call cosmlg(1)%init(); call sinmlg(1)%init()
        call da_add(pool%coslg, 0.0_DP, cosmlg(1))
        call da_add(pool%sinlg, 0.0_DP, sinmlg(1))
    
        if (gf%ncs < 2) then
            call cosmlg(1)%destroy(); call sinmlg(1)%destroy()
            call pool%destroy()
            return
        end if
        if (gf%ncs > ndeg_max) stop
        
        ! 调用 plmx_da 计算 plm 数组
        call plmx_da(gf%ncs, pool%zr, plm, pool)
        
        ! k_v = (0,0,1)
        call vec_mul(0.0_DP, pool%k_v, pool%k_v)
        call da_add(pool%k_v%elements(3), 1.0_DP, pool%k_v%elements(3))
        
        ! g_v = (-sinlg, coslg, 0)
        call da_mul(pool%sinlg, -1.0_DP, pool%g_v%elements(1))
        call da_add(pool%coslg, 0.0_DP, pool%g_v%elements(2))
        call da_mul(pool%g_v%elements(3), 0.0_DP, pool%g_v%elements(3))
        
        ! eta = sqrt(1.0_DP - zr**2)  (使用 subroutine 版本避免临时句柄泄漏)
        call da_mul(pool%zr, pool%zr, pool%tmp_da1)
        call da_mul(pool%tmp_da1, -1.0_DP, pool%tmp_da2)
        call da_add(pool%tmp_da2, 1.0_DP, pool%tmp_da3)
        call da_sqrt_sub(pool%tmp_da3, pool%eta)

        if (pool%eta%cons() < 1e-7_DP) stop
    
        do l = 2, gf%ncs
            call cosmlg(l)%init()
            call sinmlg(l)%init()
            call da_mul(cosmlg(l), 0.0_DP, cosmlg(l))
            call da_mul(sinmlg(l), 0.0_DP, sinmlg(l))

            if (l == 2) then
                ! cosmlg(l) = 2.0_DP*coslg*cosmlg(l-1) - 1.0_DP
                call da_mul(pool%coslg, cosmlg(l-1), pool%tmp_da1)
                call da_mul(pool%tmp_da1, 2.0_DP, pool%tmp_da2)
                call da_sub(pool%tmp_da2, 1.0_DP, cosmlg(l))
                ! sinmlg(l) = 2.0_DP*coslg*sinmlg(l-1)
                call da_mul(pool%coslg, sinmlg(l-1), pool%tmp_da1)
                call da_mul(pool%tmp_da1, 2.0_DP, sinmlg(l))
            else
                ! cosmlg(l) = 2.0_DP*coslg*cosmlg(l-1) - cosmlg(l-2)
                call da_mul(pool%coslg, cosmlg(l-1), pool%tmp_da1)
                call da_mul(pool%tmp_da1, 2.0_DP, pool%tmp_da2)
                call da_sub(pool%tmp_da2, cosmlg(l-2), cosmlg(l))
                ! sinmlg(l) = 2.0_DP*coslg*sinmlg(l-1) - sinmlg(l-2)
                call da_mul(pool%coslg, sinmlg(l-1), pool%tmp_da1)
                call da_mul(pool%tmp_da1, 2.0_DP, pool%tmp_da2)
                call da_sub(pool%tmp_da2, sinmlg(l-2), sinmlg(l))
            end if
            
            do m = 1, l
                ! w11 = sqrt((l+m+1.0_DP)*(l-m))  (实数标量，不需要 DA)
                ! w21 = real(m, DP)*zr/eta
                call da_mul(pool%zr, real(m, DP), pool%tmp_da1)
                call da_div(pool%tmp_da1, pool%eta, pool%w21_da)
                
                ! dplm = (w11*plm(l,m+1) - w21*plm(l,m))/eta
                ! 其中 w11 = sqrt((l+m+1)*(l-m)) 是实数
                call da_mul(plm(l,m+1), sqrt(real((l+m+1)*(l-m), DP)), pool%tmp_da1)
                call da_mul(pool%w21_da, plm(l,m), pool%tmp_da2)
                call da_sub(pool%tmp_da1, pool%tmp_da2, pool%tmp_da3)
                call da_div(pool%tmp_da3, pool%eta, pool%dplm)
                
                ! w1 = (rm/r)**l/r**3 * (((l+1.0_DP)*plm(l,m)+zr*dplm)*dr_da_tmp - r*dplm*k_v)
                ! 先计算 (rm/r)**l
                ! 注意: da_div 不支持 real/DA，先将 rm 转为 DA
                call da_mul(pool%tmp_da1, 0.0_DP, pool%tmp_da1)
                call da_add(pool%tmp_da1, rm, pool%tmp_da1)
                call da_div(pool%tmp_da1, pool%r, pool%tmp_da1)
                ! 使用 da_pow_int_sub 避免临时句柄泄漏
                call da_pow_int_sub(pool%tmp_da1, l, pool%tmp_da2)
                ! 保存 (rm/r)**l 到 ratio_tmp，供后续 w2 计算使用
                call da_add(pool%tmp_da2, 0.0_DP, pool%ratio_tmp)
                ! 再除以 r**3
                call da_mul(pool%r, pool%r, pool%tmp_da3)
                call da_mul(pool%tmp_da3, pool%r, pool%tmp_da3)
                call da_div(pool%tmp_da2, pool%tmp_da3, pool%tmp_da4)
                ! 计算 ((l+1)*plm(l,m) + zr*dplm)
                call da_mul(plm(l,m), real(l+1, DP), pool%tmp_da5)
                call da_mul(pool%zr, pool%dplm, pool%tmp_da6)
                call da_add(pool%tmp_da5, pool%tmp_da6, pool%tmp_da7)
                ! 乘以 dr_da_tmp
                call vec_mul(pool%tmp_da7, pool%dr_da_tmp, pool%tmp_vec1)
                ! 计算 r*dplm*k_v
                call da_mul(pool%r, pool%dplm, pool%tmp_da8)
                call vec_mul(pool%tmp_da8, pool%k_v, pool%tmp_vec2)
                ! 相减
                call vec_sub(pool%tmp_vec1, pool%tmp_vec2, pool%tmp_vec3)
                ! 乘以 (rm/r)**l/r**3
                call vec_mul(pool%tmp_da4, pool%tmp_vec3, pool%w1)
                
                ! w11 = clm(l,m)*cosmlg(m) + slm(l,m)*sinmlg(m)  (实数标量)
                ! 注意: w11 是实数，这里用 pool%w11_da 存储 DA 结果
                call da_mul(cosmlg(m), clm(l,m), pool%tmp_da1)
                call da_mul(sinmlg(m), slm(l,m), pool%tmp_da2)
                call da_add(pool%tmp_da1, pool%tmp_da2, pool%w11_da)
                
                ! w2 = real(m, DP)/r2/r*(rm/r)**l * plm(l,m)*g_v
                ! 计算 m/r2/r
                call da_mul(pool%w11_da, 0.0_DP, pool%tmp_da1)  ! 复用 w11_da 作为 r2
                call da_add(pool%w11_da, 0.0_DP, pool%tmp_da1)  ! 暂存 r2
                ! 实际上 r2 是 pool%w11_da 但被覆盖了，需要重新计算
        ! 重新计算 r2 = sqrt(dr_da_tmp(1)**2 + dr_da_tmp(2)**2)
        call da_mul(pool%dr_da_tmp%elements(1), pool%dr_da_tmp%elements(1), pool%tmp_da1)
        call da_mul(pool%dr_da_tmp%elements(2), pool%dr_da_tmp%elements(2), pool%tmp_da3)
        call da_add(pool%tmp_da1, pool%tmp_da3, pool%tmp_da4)
        call da_sqrt_sub(pool%tmp_da4, pool%tmp_da5)
        ! 现在 tmp_da5 = r2
        call da_mul(pool%tmp_da5, pool%r, pool%tmp_da6)
        ! da_div 不支持 real/DA，先将 real(m) 转为 DA
        call da_mul(pool%tmp_da7, 0.0_DP, pool%tmp_da7)
        call da_add(pool%tmp_da7, real(m, DP), pool%tmp_da7)
        call da_div(pool%tmp_da7, pool%tmp_da6, pool%tmp_da7)
        ! 乘以 (rm/r)**l (保存在 pool%ratio_tmp 中)
        call da_mul(pool%tmp_da7, pool%ratio_tmp, pool%tmp_da8)
        ! 乘以 plm(l,m)
        call da_mul(pool%tmp_da8, plm(l,m), pool%tmp_da1)
        ! 乘以 g_v
        call vec_mul(pool%tmp_da1, pool%g_v, pool%w2)
                
                ! w21 = clm(l,m)*sinmlg(m) - slm(l,m)*cosmlg(m)  (实数标量)
                call da_mul(sinmlg(m), clm(l,m), pool%tmp_da1)
                call da_mul(cosmlg(m), slm(l,m), pool%tmp_da2)
                call da_sub(pool%tmp_da1, pool%tmp_da2, pool%w21_da)
                
                ! flm = flm - (w1*w11 + w2*w21)
                ! w1*w11 (向量乘实数)
                call vec_mul(pool%w11_da, pool%w1, pool%tmp_vec4)
                ! w2*w21 (向量乘实数)
                call vec_mul(pool%w21_da, pool%w2, pool%tmp_vec5)
                ! 相加
                call vec_add(pool%tmp_vec4, pool%tmp_vec5, pool%tmp_vec6)
                ! flm = flm - tmp_vec6
                call vec_sub(flm, pool%tmp_vec6, pool%tmp_vec7)
                ! call vec_add(pool%tmp_vec7, 0.0_DP*flm, flm)
                flm = pool%tmp_vec7
            end do
        end do
        
        ! flm = flm * gmm
        call vec_mul(gmm, flm, pool%tmp_vec1)
        ! call vec_add(pool%tmp_vec1, 0.0_DP*flm, flm)
        flm =  pool%tmp_vec1

        ! 清理 cosmlg, sinmlg 数组
        do l = 1, gf%ncs
            call cosmlg(l)%destroy()
            call sinmlg(l)%destroy()
        end do
        ! 清理 plm 数组
        do l = 1, size(plm, 1)
            do m = 1, size(plm, 2)
                call plm(l,m)%destroy()
            end do
        end do
        
        ! 销毁临时变量池
        call pool%destroy()
    end subroutine f_tesseral_da
    
    ! =========================================================
    ! DA 版本：Legendre 多项式（显式运算）
    ! =========================================================
    subroutine plx_da(n, zr, pl, pool)
        integer, intent(in) :: n
        type(DA), intent(in) :: zr
        type(DA), intent(out) :: pl(n)
        type(GravityTempPool), intent(inout) :: pool
        integer :: l
        real(DP) :: w1, w2, l2
        
        do l = 1, n
            call pl(l)%init()
            call da_mul(pl(l), 0.0_DP, pl(l))
        end do

        ! pl(1) = sqrt(3.0_DP) * zr
        call da_mul(zr, sqrt(3.0_DP), pl(1))
        
        ! pl(2) = sqrt(5.0_DP)/2.0_DP * (3.0_DP*zr*zr - 1.0_DP)
        call da_mul(zr, zr, pool%tmp_da1)
        call da_mul(pool%tmp_da1, 3.0_DP, pool%tmp_da2)
        call da_sub(pool%tmp_da2, 1.0_DP, pool%tmp_da3)
        call da_mul(pool%tmp_da3, sqrt(5.0_DP)/2.0_DP, pl(2))
        
        if (n < 3) return
        do l = 3, n
            l2 = 2.0_DP * real(l, DP)
            w1 = sqrt((l2+1.0_DP)/(l2-1.0_DP))
            w2 = sqrt((l2-1.0_DP)/(l2-3.0_DP))
            
            ! pl(l) = w1 * ((2.0_DP-1.0_DP/real(l, DP))*zr*pl(l-1) - w2*(1.0_DP-1.0_DP/real(l, DP))*pl(l-2))
            ! 计算 (2.0-1.0/l)*zr*pl(l-1)
            call da_mul(zr, pl(l-1), pool%tmp_da1)
            call da_mul(pool%tmp_da1, (2.0_DP-1.0_DP/real(l, DP)), pool%tmp_da2)
            ! 计算 w2*(1.0-1.0/l)*pl(l-2)
            call da_mul(pl(l-2), (1.0_DP-1.0_DP/real(l, DP)), pool%tmp_da3)
            call da_mul(pool%tmp_da3, w2, pool%tmp_da4)
            ! 相减
            call da_sub(pool%tmp_da2, pool%tmp_da4, pool%tmp_da5)
            ! 乘以 w1
            call da_mul(pool%tmp_da5, w1, pl(l))
        end do
    end subroutine plx_da

    ! =========================================================
    ! DA 版本：连带 Legendre 多项式（显式运算）
    ! =========================================================
    subroutine PlmX_da(n, zr, plm, pool)
        integer, intent(in) :: n
        type(DA), intent(in) :: zr
        type(DA), intent(out) :: plm(n,n+1)
        type(GravityTempPool), intent(inout) :: pool
        integer :: l, m
        real(DP) :: w1, w2, l2
        
        ! eta = sqrt(1.0_DP - zr*zr)  (使用 subroutine 版本避免临时句柄泄漏)
        call da_mul(zr, zr, pool%tmp_da1)
        call da_mul(pool%tmp_da1, -1.0_DP, pool%tmp_da2)
        call da_add(pool%tmp_da2, 1.0_DP, pool%tmp_da3)
        call da_sqrt_sub(pool%tmp_da3, pool%eta)

        ! 初始化 plm 数组
        do m = 1, n+1
            do l = 1, n
                call plm(l,m)%init()
                call da_mul(plm(l,m), 0.0_DP, plm(l,m))
            end do
        end do
        
        ! plm(1,1) = sqrt(3.0_DP) * eta
        call da_mul(pool%eta, sqrt(3.0_DP), plm(1,1))
        
        if (n < 2) return

        ! plm(2,1) = sqrt(5.0_DP) * zr * plm(1,1)
        call da_mul(zr, plm(1,1), pool%tmp_da1)
        call da_mul(pool%tmp_da1, sqrt(5.0_DP), plm(2,1))
        
        ! plm(2,2) = sqrt(5.0_DP)/2.0_DP * eta * plm(1,1)
        call da_mul(pool%eta, plm(1,1), pool%tmp_da1)
        call da_mul(pool%tmp_da1, sqrt(5.0_DP)/2.0_DP, plm(2,2))
        
        if (n < 3) return

        do l = 3, n
            l2 = 2.0_DP * real(l, DP)
            do m = 1, l-1
                w1 = (l2+1.0_DP)*(l2-1.0_DP) / (real(l+m, DP)*real(l-m, DP))
                w2 = (l2+1.0_DP)*(real(l-1+m, DP)*real(l-1-m, DP)) / ((l2-3.0_DP)*real(l+m, DP)*real(l-m, DP))
                
                ! plm(l,m) = sqrt(w1)*zr*plm(l-1,m) - sqrt(w2)*plm(l-2,m)
                call da_mul(zr, plm(l-1,m), pool%tmp_da1)
                call da_mul(pool%tmp_da1, sqrt(w1), pool%tmp_da2)
                call da_mul(plm(l-2,m), sqrt(w2), pool%tmp_da3)
                call da_sub(pool%tmp_da2, pool%tmp_da3, plm(l,m))
            end do
            ! plm(l,l) = sqrt((l2+1.0_DP)/l2) * eta * plm(l-1,l-1)
            call da_mul(pool%eta, plm(l-1,l-1), pool%tmp_da1)
            call da_mul(pool%tmp_da1, sqrt((l2+1.0_DP)/l2), plm(l,l))
        end do
    end subroutine PlmX_da

end module pod_gravity_model_module
