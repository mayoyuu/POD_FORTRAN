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
    
contains

    subroutine read_gravity_field(gf)
        class(gravity_field), intent(in) :: gf
        character(100) :: filename
        integer :: i1, i2, nstat, ios
        real(DP), dimension(2) :: tmp
        
        ! 修改为跨平台的相对路径
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
   
    subroutine f_tesseral(gf, flm)
        class(gravity_field), intent(in) :: gf
        integer :: l, m
        real(DP), dimension(gf%ncs) :: cosmlg, sinmlg
        real(DP), dimension(gf%ncs, gf%ncs+1) :: plm
        real(DP) :: coslg, sinlg, r, r2, zr, w11, w21, eta, dplm
        real(DP), dimension(3) :: k_v, g_v, w1, w2, flm
        
        ! --- 安全检查提前 ---
        if (gf%ncs < 1) return ! 如果没有阶数，直接返回
        
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

    !!!! DA版本的引力加速度计算接口
    subroutine f_zonal_da(gf, fl)
        class(gravity_field), intent(in) :: gf
        type(AlgebraicVector), intent(out) :: fl
        type(AlgebraicVector) :: h, dr_da_tmp
        integer :: l
        type(DA), dimension(gf%ncs) :: pl
        type(DA) :: r, r3, rl3, rg2, zr, tmp2, tmp3, ratio_tmp
        real(DP) :: tmp1
        
        if (gf%cen_body == 3) then
            cl0 = cl0e; rm = r_e; gmm = g_e
        end if
        if (gf%cen_body == 10) then
            cl0 = cl0m; rm = r_m; gmm = g_m
        end if
        
        call fl%init(3)
        call h%init(3)

        fl  = 0.0_DP

        dr_da_tmp = gf%dr_da
        r = dr_da_tmp%norm2()

        zr = dr_da_tmp%elements(3) / r
        r3 = r**3
        rg2 = dr_da_tmp%elements(1)**2 + dr_da_tmp%elements(2)**2

        h%elements(1) = dr_da_tmp%elements(1)*dr_da_tmp%elements(3)
        h%elements(2) = dr_da_tmp%elements(2)*dr_da_tmp%elements(3)
        h%elements(3) = -rg2

        if (gf%ncs < 2) return
        call plx_da(gf%ncs, zr, pl)
        
        do l = 2, gf%ncs
            rl3 = ((r/rm)**l) * r3
            tmp1 = sqrt((2.0_DP*l+1.0_DP)/(2.0_DP*real(l, DP)-1.0_DP))
            tmp2 = (l+1.0_DP) * pl(l)
            
            ratio_tmp = rg2/r
            if (ratio_tmp%cons() > 1.0e-14_DP) then 
                tmp3 = real(l, DP)*r/rg2 * (tmp1*pl(l-1) - zr*pl(l))
            else 
                tmp3 = 1.0_DP
            end if
            fl = fl - cl0(l) * (tmp2*dr_da_tmp + tmp3*h) / rl3
        end do
        fl = fl * gmm 
    end subroutine f_zonal_da

    subroutine f_tesseral_da(gf, flm)
        class(gravity_field), intent(in) :: gf
        type(AlgebraicVector) :: dr_da_tmp
        type(AlgebraicVector), intent(out) :: flm
        

        integer :: l, m
        type(DA), dimension(gf%ncs) :: cosmlg, sinmlg
        type(DA), dimension(gf%ncs, gf%ncs+1) :: plm
        type(DA) :: coslg, sinlg, r, r2, zr, w11, w21, eta, dplm
        type(AlgebraicVector) :: k_v, g_v, w1, w2
        
        ! --- 安全检查提前 ---
        if (gf%ncs < 1) return ! 如果没有阶数，直接返回
        
        if (gf%cen_body == 3) then
            clm = clme; slm = slme
        else 
            clm = clmm; slm = slmm
        end if
        
        call flm%init(3)
        flm = 0.0_DP

        dr_da_tmp = gf%dr_da
        r = dr_da_tmp%norm2()
        r2 = sqrt(dr_da_tmp%elements(1)**2 + dr_da_tmp%elements(2)**2)
        zr = dr_da_tmp%elements(3) / r
        coslg = dr_da_tmp%elements(1) / r2
        sinlg = dr_da_tmp%elements(2) / r2

        call cosmlg(1)%init(); call sinmlg(1)%init()
        cosmlg(1) = coslg; sinmlg(1) = sinlg
    
        if (gf%ncs < 2) return
        if (gf%ncs > ndeg_max) stop
        
        call plmx_da(gf%ncs, zr, plm)
        
        call k_v%init(3)
        call g_v%init(3)
        k_v = 0.0_DP; k_v%elements(3) = 1.0_DP
        g_v%elements(1) = -sinlg; g_v%elements(2) = coslg; g_v%elements(3) = 0.0_DP
        eta = sqrt(1.0_DP - zr**2)

        if (eta%cons() < 1e-7_DP) stop
    
        do l = 2, gf%ncs
            call cosmlg(l)%init()
            call sinmlg(l)%init()
            cosmlg(l) = 0.0_DP; sinmlg(l) = 0.0_DP

            if (l == 2) then
                cosmlg(l) = 2.0_DP*coslg*cosmlg(l-1) - 1.0_DP
                sinmlg(l) = 2.0_DP*coslg*sinmlg(l-1)
            else
                cosmlg(l) = 2.0_DP*coslg*cosmlg(l-1) - cosmlg(l-2)
                sinmlg(l) = 2.0_DP*coslg*sinmlg(l-1) - sinmlg(l-2)
            end if
            
            do m = 1, l
                w11 = sqrt((l+m+1.0_DP)*(l-m))
                w21 = real(m, DP)*zr/eta
                dplm = (w11*plm(l,m+1) - w21*plm(l,m))/eta
                w1 = (rm/r)**l/r**3 * (((l+1.0_DP)*plm(l,m)+zr*dplm)*dr_da_tmp - r*dplm*k_v)
                w11 = clm(l,m)*cosmlg(m) + slm(l,m)*sinmlg(m)
                w2 = real(m, DP)/r2/r*(rm/r)**l * plm(l,m)*g_v
                w21 = clm(l,m)*sinmlg(m) - slm(l,m)*cosmlg(m)
                flm = flm - (w1*w11 + w2*w21)
            end do
        end do
        flm = flm * gmm
    end subroutine f_tesseral_da
    
    subroutine plx_da(n, zr, pl)
        integer, intent(in) :: n
        type(DA), intent(in) :: zr
        type(DA), intent(out) :: pl(n)
        integer :: l
        real(DP) :: w1, w2, l2
        
        do l = 1, n
            call pl(l)%init()
            pl(l) = 0.0_DP
        end do

        pl(1) = sqrt(3.0_DP) * zr
        pl(2) = sqrt(5.0_DP)/2.0_DP * (3.0_DP*zr*zr - 1.0_DP)
        if (n < 3) return
        do l = 3, n
            l2 = 2.0_DP * real(l, DP)
            w1 = sqrt((l2+1.0_DP)/(l2-1.0_DP))
            w2 = sqrt((l2-1.0_DP)/(l2-3.0_DP))
            pl(l) = w1 * ((2.0_DP-1.0_DP/real(l, DP))*zr*pl(l-1) - w2*(1.0_DP-1.0_DP/real(l, DP))*pl(l-2))
        end do
    end subroutine plx_da

    subroutine PlmX_da(n, zr, plm)
        integer, intent(in) :: n
        type(DA), intent(in) :: zr
        type(DA), intent(out) :: plm(n,n)
        integer :: l, m
        type(DA) :: eta
        real(DP) :: w1, w2, l2
        
        eta = sqrt(1.0_DP - zr*zr)

        do l = 1, n
            do m = 1, l
                call plm(l,m)%init()
                plm(l,m) = 0.0_DP
            end do
        end do
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
    end subroutine PlmX_da


end module pod_gravity_model_module
