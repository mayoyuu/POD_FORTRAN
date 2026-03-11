! ======================================================================

!   Gravity field of a centeral body
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  02-07-2024 (created)
!   ref: Liu Lin 《卫星轨道力学与应用》; code from Hu Shoucun
!
!   Noticed that Here defined a global varaible "gm" that may lead to repetitive naming with DE440 GM
!
! ======================================================================
!
!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================

Module gravity_module
    
    use kind_parameter, only: dp
    
    implicit none
    
    real(dp)                                :: g_m,r_m,g_e,r_e,gmm,rm
    integer, parameter                        :: ndeg_max = 100
    real(dp),dimension(ndeg_max)            :: cl0
    real(dp),dimension(100),public            :: cl0e
    real(dp),dimension(100),public            :: cl0m
    real(dp),dimension(ndeg_max,ndeg_max)    :: clm,slm
    real(dp),dimension(ndeg_max,ndeg_max)    :: clme,clmm,slme,slmm
    
    type,public:: gravity_field
        
        real(dp),dimension(3)                :: dr
        integer                                :: ncs,cen_body
        
    contains
        procedure,public,pass(gf):: read_gravity_field
        procedure,public,pass(gf):: f_zonal
        procedure,public,pass(gf):: f_tesseral
        
    end type gravity_field
    
    contains
    
! ======================================================================

!   Read the gravity field of a centeral body
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  12-07-2024 (created)
!
! ======================================================================
!
!     Input files:     Gravity filed model: Gm05C.GEO, jggrx_0900d_sha.tab
!                        
!     Output files:    
!
! ======================================================================
!
!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================

    subroutine read_gravity_field (gf)
    
    
        implicit none
        
        class(gravity_field),intent(in)        :: gf
        character(100)                                :: filename
        integer                                        :: ip,i1,i2,nstat,k,ios
        real(dp),dimension(2)                        :: tmp
        if (gf%cen_body == 3)  filename = 'gravity_model\GGM05C.GEO'
        if (gf%cen_body == 10) filename = 'gravity_model\gggrx_0660pm_sha.tab'     
        cl0 = 0.0_dp
        clm = 0.0_dp
        slm = 0.0_dp
        open (101, file = filename,iostat = ios)
        ! print*,'OPEN',ios
        if(ios==0)then
        if (gf%cen_body == 3)   then
            read (101, *) g_e,r_e
            ! print*,'e',g_e,r_e
            g_e = g_e /1e9_dp
            r_e = r_e /1e3_dp
        end if
        if (gf%cen_body == 10)  then
            read (101, *) g_m,r_m
            ! print*,'m',g_m,r_m
            g_m = g_m /1e9_dp
            r_m = r_m /1e3_dp
        end if
        else
            print*,'file open failed'
        endif
        
        
        do while (.true.)
        
            read(101,*,iostat=nstat) i1,i2,tmp
            
            if (nstat /=0 ) exit
            if (i1 > gf%ncs .or. i2 > gf%ncs ) exit
            if (i1 > ndeg_max .or. i2 > ndeg_max ) exit
            
            if (i2 == 0) then
            
                !cl0(i1) = tmp(1)
                
                if (gf%cen_body == 3) cl0e(i1) = tmp(1)
                if (gf%cen_body == 10) cl0m(i1) = tmp(1)
            else 
            
                !clm(i1,i2) = tmp(1)
                !slm(i1,i2) = tmp(2)
                if (gf%cen_body == 3) then
                    clme(i1,i2) = tmp(1)
                    slme(i1,i2) = tmp(2)
                else
                    clmm(i1,i2) = tmp(1)
                    slmm(i1,i2) = tmp(2)
                end if
                
            end if
            
        end do    
        
        ! print*,cl0e(2),cl0m(2)
        
        
        close (101)
        
        return
        end subroutine
        


! ======================================================================
!
!   Zonal gravity field of a centeral body
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  12-07-2024 (created)
!   ref: Liu Lin 《卫星轨道力学与应用》; code from Hu Shoucun
!
! ======================================================================

!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================
    
    
    subroutine f_zonal(gf,fl)

        
        implicit none
        
        class(gravity_field),intent(in)        :: gf
        integer(4)                            :: l,m
        real(dp),dimension(gf%ncs)            :: pl
        real(dp)                            :: r,r3,rl3,zr,rg2
        real(dp)                            :: tmp1,tmp2,tmp3
        real(dp),dimension(3)                :: h,fl
        ! print*,'begin'
        if (gf%cen_body == 3) then
            cl0 = cl0e
            rm = r_e
            gmm = g_e
            ! print*,'gmm',gmm
        end if
        
        if (gf%cen_body == 10) then
            cl0 = cl0m
            rm = r_m
            gmm = g_m
        end if
        
        
        fl     = 0.0_dp
        r     = norm2(gf%dr)
        ! print*,'fl,r',fl,r
        zr     = gf%dr(3)/r
        r3     = r**3_dp
        rg2 = gf%dr(1)**2_dp+gf%dr(2)**2_dp
        h     = (/gf%dr(1)*gf%dr(3), gf%dr(2)*gf%dr(3), -rg2/)


        
        if (gf%ncs < 2) return
        
        call plx(gf%ncs, zr, pl)
        
        do l = 2, gf%ncs
            rl3 = (r/rm)**real(l)*r3
            tmp1 = ((2.0_dp*l+1_dp)/(2_dp*real(l)-1_dp))**(0.5_dp)
            tmp2 = (l+1_dp)*pl(l)
            
            if (rg2/r > 1.0D-14) then  !无奇点
                tmp3 = real(l)*r/rg2*(tmp1*pl(l-1)-zr*pl(l))
            else                       !奇点近似处理
                tmp3 = 1_dp
            end if
            
            fl = fl-cl0(l)*(tmp2*gf%dr+tmp3*h)/rl3
        end do

        fl = fl*gmm !*(rm/dsqrt(rm**3/gm))
   return
   end subroutine
   
   
! ======================================================================
!
!   Tesseral gravity field of a centeral body
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  12-07-2024 (created)
!   ref: Liu Lin 《卫星轨道力学与应用》; code from Hu Shoucun
!
! ======================================================================

!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================
            
            
        
    subroutine f_tesseral(gf,flm)
    
    
        implicit none
        
        class(gravity_field),intent(in)        :: gf
        integer(4) :: l,m
        real(dp),dimension(gf%ncs)            :: cosmlg,sinmlg
        real(dp),dimension(gf%ncs,gf%ncs+1)    :: plm
        real(dp)                            :: coslg,sinlg,r,r2,zr,w11,w21,eta,dplm
        real(dp),dimension(3)                :: k_v,g_v,w1,w2,flm
        
        if (gf%cen_body == 3) then
            clm = clme
            slm = slme
        else 
            clm = clmm
            slm = slmm
        end if
            
    
        flm = 0.0_dp
        r     = norm2(gf%dr)
        r2     = sqrt(gf%dr(1)**2_dp+gf%dr(2)**2_dp)
        zr     = gf%dr(3)/r
        coslg = gf%dr(1)/r2
        sinlg = gf%dr(2)/r2
        cosmlg(1) = coslg
        sinmlg(1) = sinlg
    
        if (gf%ncs < 2) return
        if (gf%ncs > ndeg_max) stop
        
        call plmx(gf%ncs,zr,plm)
        
        k_v = (/0_dp,0_dp,1_dp/)
        g_v(1) = -sinlg
        g_v(2) = coslg
        g_v(3) = 0_dp
        eta = sqrt(1-zr**2d0)

        if (eta < 1d-7) stop

    
        do l = 2,gf%ncs
            if (l == 2) then
                cosmlg(l) = 2d0*coslg*cosmlg(l-1)-1d0
                sinmlg(l) = 2d0*coslg*sinmlg(l-1)
                goto 10
            end if
            cosmlg(l) = 2d0*coslg*cosmlg(l-1)-cosmlg(l-2)
            sinmlg(l) = 2d0*coslg*sinmlg(l-1)-sinmlg(l-2)
10          do m = 1,l
                w11 = sqrt((l+m+1d0)*(l-m))
                w21 = m*zr/eta
                dplm = (w11*plm(l,m+1)-w21*plm(l,m))/eta
                w1 = (rm/r)**l/r**3 * (((l+1d0)*plm(l,m)+zr*dplm)*gf%dr-r*dplm*k_v)
                w11 = clm(l,m)*cosmlg(m)+slm(l,m)*sinmlg(m)
                w2 = m/r2/r*(rm/r)**l * plm(l,m)*g_v
                w21 = clm(l,m)*sinmlg(m)-slm(l,m)*cosmlg(m)
                flm = flm-(w1*w11+w2*w21)
            end do
        end do
        flm = flm*gmm
        return
    end subroutine    
        
    
    
    
    
    
    
! ======================================================================
!
!   Calculate Pl
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  12-07-2024 (created)
!   ref: Liu Lin 《卫星轨道力学与应用》; code from Hu Shoucun
!
! ======================================================================

!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================
    
    
    
    
    
    subroutine plx(n,zr,pl)
        implicit none
        integer(4) :: n,l,l2
        real(dp) :: zr,w1,w2,pl(n)
        pl = 0d0
        pl(1) = sqrt(3d0)*zr
        pl(2) = sqrt(5d0)/2d0*(3d0*zr*zr-1d0)
        if (n < 3) return
        do l = 3,n
        l2 = 2d0*l
        w1 = sqrt((l2+1d0)/(l2-1d0))
        w2 = sqrt((l2-1d0)/(l2-3d0))
        pl(l) = w1*((2d0-1d0/l)*zr*pl(l-1)-w2*(1d0-1d0/l)*pl(l-2))
    end do
    return
    end subroutine
    

! ======================================================================
!
!   Caculate Plm
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  12-07-2024 (created)
!   ref: Liu Lin 《卫星轨道力学与应用》; code from Hu Shoucun
!
! ======================================================================

!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================

    subroutine PlmX(n,zr,plm)
    implicit none
    integer(4) :: n,l,m
    real(8) :: zr,plm(n,n),eta,w1,w2,l2
    real(8),external :: factorial
    eta = sqrt(1d0-zr**2d0)
    plm = 0d0
    plm(1,1) = sqrt(3d0)*eta
    plm(2,1) = sqrt(5d0)*zr*plm(1,1)
    plm(2,2) = sqrt(5d0)/2d0*eta*plm(1,1)
    if (n < 3) return
    do l = 3,n
        l2 = 2d0*l
        do m = 1,l-1
            w1 = (l2+1d0)*(l2-1d0)/(l+m)/(l-m)
            w2 = (l2+1d0)*(l-1d0+m)*(l-1d0-m)/(l2-3d0)/(l+m)/(l-m)
            plm(l,m) = sqrt(w1)*zr*plm(l-1,m)-sqrt(w2)*plm(l-2,m)
        end do
        plm(l,l) = sqrt((l2+1d0)/l2)*eta*plm(l-1,l-1)
    end do
    return
    end subroutine

    
    
    
    
end module
