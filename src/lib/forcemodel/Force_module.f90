! ======================================================================

!   Forces acting on the object.
!
! ======================================================================
!
!   Author: Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!
!   Update history:  02-07-2024 (created)
!   ref: Liu Lin 《卫星轨道力学与应用》; code from Hu Shoucun
!
!
! ======================================================================
!
!     Variable  I/O  Description
!     --------  ---  --------------------------------------------------
!
! ======================================================================

    Module force_module
    
        use kind_parameter, only: dp
        use gravity_module
        USE M_DE440
        use constants_module
        implicit none
        
        public::gravity_const
        public::yhc
        public::grav_planets
        
        real(dp),dimension(11)                    :: gm_cons
        real(dp)                                :: n_nsp
        real(dp)                                :: tdb0,tdb




    contains
    
    subroutine gravity_const
    
        implicit none
        
        ! integer(4)                                :: i
        integer(4)                                :: drm!,dm,ip
        ! REAL(DP)                                    ::AU,GM(11)
        
        drm = 5
        
        ! do i = 1, 11
        

        !     if (i == 10)  then
        !         ip = 301
        !         call bodvcd(ip,'GM',drm,dm,gm_cons(i))
                
        !     else if (i == 11) then
        !         ip = 10
        !         call bodvcd(ip,'GM',drm,dm,gm_cons(i))
        !     else 
        !         ip = i
        !         call bodvcd(ip,'GM',drm,dm,gm_cons(i))
        !     end if
            
        !     ! print*,i, gm_cons(i)
            
            
        ! end do
        ! CALL CONSTANTS(AU,GM)
        gm_cons = GM    
    end subroutine
    
!!=============================================================================
!! 	accoording to orbitprop-main's rkf78_module; t1-real-time;x-the state vector;y-the output F
    subroutine yhc (t1,x,y)

        use gravity_module
        
        implicit none
        real(dp),dimension(6)                    :: x,y,rp,rxp
        real(dp)                                :: r,drp,drxp,lt,t1
        integer, dimension(11)                     :: i_planets
        integer                                    :: i,ndeg,n
        real(dp),dimension(11,3)                :: fp
        
        y = 0.0_dp
        
        y(1:3) = x(4:6)
        
        fp = 0.0_dp
        
        i_planets = 0
        i_planets(3) = 1
        i_planets(10) = 1
        i_planets(11) = 1
        
        do i = 1,11
        
            if (i_planets(i) > 0) then
            
                ndeg = 0
            
                if (i == 3 .or. i == 10) ndeg = int(n_nsp)
                !if (i == 3) ndeg = 10

                call grav_planets(t1, x, i,fp,ndeg)

                !print*,i,fp(i,1:3)
                
                y(4:6) = y(4:6) + fp(i,1:3)
            
            end if
            
        end do

        
        
    end subroutine
    
    
    subroutine yhc1 (t1,x,y,n)
        
        implicit none
        real(dp),dimension(6)                    :: x,y,rp,rxp
        real(dp)                                :: r,drp,drxp,lt,t1,t
        integer, dimension(11)                     :: i_planets
        integer                                    :: i,ndeg,n
        real(dp),dimension(11,3)                :: fp
        
        y = 0.0_dp
        
        y(1:3) = x(4:6)
        
        fp = 0.0_dp
        
        i_planets = 0
        i_planets(3) = 1
        i_planets(10) = 1
        i_planets(11) = 1
        
        do i = 1,11
        
            if (i_planets(i) > 0) then
            
                ndeg = 0
            
                if (i == 3 .or. i == 10) ndeg = 100
            
                call grav_planets(t1, x, i,fp,ndeg)
                
                ! print*,i,fp(i,1:3)
                
                y(4:6) = y(4:6) + fp(i,1:3)
            
            end if
            
            
        end do

        
        
    end subroutine
        
    
    subroutine grav_planets(t1, x, i,fp,ndeg)
    
        implicit none
        
    
        type(gravity_field)                        :: fsce,fscm
        real(dp),dimension(6),intent(in)        :: x
        integer,intent(in)                        :: i,ndeg
        integer                                    :: ib,np
        real(dp),dimension(3)                    :: rp,rxp
        real(dp)                                :: r,drp,drxp,lt,t1
        real(dp),dimension(3)                    :: fl,flm,xb,fl0,flm0
        real(dp),dimension(11,3)                :: fp
        real(dp),dimension(3,3)                    :: xform
        character(4)                            :: ibody,timstr
        logical                                    :: found
        
        r = norm2(x(1:3))        

        if (i ==  3)    then
                
            fp(i,1:3) = -gm_cons(i)*x(1:3)/r**3
            !fp(i,1:3) = -398600.44150_dp*x(1:3)/r**3
            
            if (ndeg > 0) then 
            
                fsce%ncs = ndeg
                fsce%cen_body = i

                xform = 0.0_dp
                call pxform ( 'J2000', 'IAU_EARTH', t1+tdb0, xform)
                call mxv(xform, x(1:3), xb)
                fsce%dr = xb(1:3)
                
                ! call fsce%read_gravity_field()
                call fsce%f_zonal(fl)
                
                call fsce%f_tesseral(flm)
                
                xform = 0.0_dp
                call pxform (  'IAU_EARTH', 'J2000',t1+tdb0, xform)
                call mxv(xform, fl, fl0)
                call mxv(xform, flm, flm0)
                
                                
                fp(i,:) = fp(i,:) + fl0 + flm0
                
                
                
            end if
            
        else 
            
            if (i == 10) then
                ibody = 'MOON'
                gm_cons(i) = 4902.800218526379882812_dp
            else if (i == 11) then
                ibody = 'SUN'
            else
                ib = i*100+99
                call bodc2n(ib, ibody, found)
            end if


            call spkpos(ibody, t1+tdb0, 'j2000', 'none', 'earth', rp, lt)
            
            drp = norm2(rp(1:3))
            
            rxp(1:3) = x(1:3)-rp(1:3)
            
            drxp = norm2(rxp)
            
            fp(i,1:3) = -gm_cons(i)*(rxp(1:3)/drxp**3+rp(1:3)/drp**3)
            
            if (i == 10 .and. ndeg > 0) then 
                        
            
                fscm%ncs = ndeg
                fscm%cen_body = i
                xform = 0.0_dp
            
                call pxform ( 'J2000', 'MOON_PA', t1+tdb0, xform)
                call mxv(xform, rxp, xb)
                fscm%dr = xb
            
                call fscm%f_zonal(fl)
                call fscm%f_tesseral(flm)
            
                !xform = 0.0_dp
                call pxform ( 'MOON_PA', 'J2000', t1+tdb0, xform)
            
                call mxv(xform, fl, fl0)
                call mxv(xform, flm, flm0)
            
                fp(i,:) = fp(i,:) + fl0 + flm0  
            
            
            
            !test lunar NSG 
            !do np = 1,100
!
            !    fscm%ncs = np !ndeg
            !    fscm%cen_body = i
            !    xform = 0.0_dp
            !    
            !    call pxform ( 'J2000', 'MOON_PA', t1+tdb0, xform)
            !    call mxv(xform, rxp, xb)
!
            !    !xb = 0.0_dp  !!!
            !    !xb(1) = 1000_dp !!!
            !    !xb(2) = 900_dp
            !    !xb(3) = 800_dp
!
            !    if (np == 0) then
            !        call pxform ( 'MOON_PA','J2000', t1+tdb0, xform)
            !        call mxv(xform, xb, rxp)
            !        print*,rxp+rp
            !    end if
!
            !    fscm%dr = xb
            !    
            !    !call fscm%read_gravity_field()
            !    call fscm%f_zonal(fl)
            !    call fscm%f_tesseral(flm)
            !    
            !    xform = 0.0_dp
!
            !    call pxform ( 'MOON_PA', 'J2000', t1+tdb0, xform)
!
            !    if (np == 1) then
            !        print*, "TM"
            !        print*, xform
            !    end if
            !    
            !    call mxv(xform, fl+flm, fl0)
            !    !call mxv(xform, flm, flm0)
            !    
            !    print*,fp(10,1:3),fl0,flm0        !    
            !    fp(i,:) = fp(i,:) + fl0 + flm0
            !    print*, "Lunar gravity field:", np
            !    print*, "LPAF:",(fl+flm)*1e6_dp
            !    print*, "GCRF:",(fl0)*1e6_dp !(fl0+flm0)*1e6_dp
!
            !end do
            
                
            end if
        end if
        
        
        
            
        return
    end subroutine

end module
