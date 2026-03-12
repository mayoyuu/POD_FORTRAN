!!========================================================================
!! the propagate orbit module in EPH,considering EMS and HFEM     
!! auther:shi zhipeng
!!
!! the propagate is in GCRS,and the units is km and km/s
!! using the Spice
!!
!! update history: 2024.6.21 created
!!                 2024.11.13 add the HFEM
!!                 2024.11.14 add the condition of propagate time<0
!!                 
!! ref : 《航天器定轨理论与应用》
!!    
module propagateOrbit_EPH_module
        use kind_parameter ,only:dp
        use M_DE440
        use math_module
        use force_module
        use tools_module
        use transform_coordinates_module        !! only for propagate_ephemeris_segments 
        use constants_module

        implicit none

        public

        public computeF_J2000                       !! compute F
        public computeDF_J2000                      !! compute DF
        public rkf78                                !! compute 3rd order analytical solution halo orbits in crtbp
        public propagateOrbit_EPH                   !! the integrator of rkf78
        public propagateOrbit_EPH_fixedstep
        ! public propagate_ephemeris_segments         !! 输出星历
       
        
        contains

!!========================================================================
!! the right function of J2000 equation of motion: 
!! auther:shi zhipeng
!!
!! update history: 2024.6.21 created
!!                 2024.11.18 add the HFEM part
!! ref : 《航天器动力学引论》
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system
!! tMJD       I   MJD时间
!! model      I   选取的模型
!! sv         I   the state vector 
!! dsv        O   the differential vector of state vector
!!
    subroutine computeF_J2000(model,tMJD,sv,dsv)
        implicit none

        real(dp),intent(in)::tMJD
        character(len=*),intent(in)::model
        real(dp),dimension(6),intent(in):: sv
        real(dp),dimension(6),intent(out)::dsv

        !! parameter
        !type(gravity_field) ::gfs
        real(dp)::unitim,ET,DSVasp(6)
        integer::i
        real(dp)::SVM(6),SVS(6),dist,F0(3),F2(3),POS(3)
        real(dp)::SVDE(6,11)
        real(dp)::PNE(6),PG(6),SRP(6)

        forall(i=1:3) dsv(i) = sv(i+3)
        POS = sv(1:3)
        
        ! call gravity_const
        ! gfs%ncs = 100 ! order and degree to be loaded
        ! gfs%cen_body = 10
        ! call gfs%read_gravity_field()
        ! gfs%cen_body = 3
        ! call gfs%read_gravity_field()
        ! n_nsp = 30  ! order and degree considred in the model(10/20/30 doesnt matter)'

        SELECT CASE (model)
        case ('EMS')
            call DE440(tMJD,10,3,SVM)
            call DE440(tMJD,11,3,SVS)
            dist = norm2(POS)
            F0 = -GM(3)*POS/dist**3
            F2 = -GM(10)*((POS-SVM(1:3))/norm2((POS-SVM(1:3)))**3+SVM(1:3)/norm2(SVM(1:3))**3)+&        !! change for GM
                 (-GM(11)*((POS-SVS(1:3))/norm2((POS-SVS(1:3)))**3+SVS(1:3)/norm2(SVS(1:3))**3))

            dsv(4:6) = F0+F2
        case ('HFEM')
            SVDE=0.0_dp
            ET = unitim(tMJD+2400000.5_dp,'JDTDB','ET')
            do i =1,11
                call DE440(tMJD,nameindex_440(i),nameindex_440(3),SVDE(:,i))
            enddo
            
            !! 大行星引力
            call f_PG(SVDE,sv,PG)
            
            !! 地月非球型和太阳第三体引力
            call yhc(ET,sv,DSVasp)
            
            !! 太阳光压
            call f_SRP(SVDE,SV,SRP)
            
            !! 相对论摄动
            call f_PNE(SVDE,sv,PNE)

            do i =1,3
                ! print*,'右函数对比',sum(Fall(i,:)) 
                dsv(i+3) = DSVasp(i+3)+PG(i+3)+SRP(i+3)+PNE(i+3) 
            enddo
        case default
            print*, "Unknown model: ", model
        end SELECT
        ! print*,'dsv',dsv
        ! print*,'SRP',SRP(4:6)
        ! print*,'PNE',PNE(4:6)
        ! stop
        end subroutine computeF_J2000
!!--------------------------------------------------------------------------------------------
!! 计算个大行星引力，日地月由yhc计算
!!
        subroutine f_PG(SVDE,sv,dsv)
            implicit none
            REAL(DP),intent(in)     ::SVDE(6,11),sv(6)
            real(dp),intent(out)    ::dsv(6)

            integer                 :: i
            real(dp)                :: POS(3),VEL(3)
            real(dp)                :: Fall(3,11)
            POS = sv(1:3)
            VEL = sv(4:6)
            do i =1,11
                if(i==3 .or. i==10 .or. i==11)then
                    Fall(:,i) = 0.0_dp
                    ! print*,i
                else            
                    Fall(:,i) = -GM(i)*((POS-SVDE(1:3,i))/norm2((POS-SVDE(1:3,i)))**3&
                    +SVDE(1:3,i)/norm2(SVDE(1:3,i))**3)
                endif
            enddo
            dsv(1:3) = VEL
            do i=1,3
                dsv(i+3) = sum(Fall(i,:))
            enddo

        end subroutine f_PG
!!-----------------------------------------------------------------------------------------------
!! 计算光压效应
!!
        subroutine f_SRP(SVDE,sv,dsv,Cr,SMR,RP)
            implicit none
            REAL(DP),intent(in)     ::SVDE(6,11),sv(6)
            real(dp),intent(out)    ::dsv(6)
            REAL(DP),optional,intent(in)        ::Cr,SMR,RP
            REAL(DP)                ::Cr0,SMR0,RP0,SVS(6),SRP(3)      !! 默认的光压常数
            real(dp)                :: POS(3),VEL(3)
            POS = sv(1:3)
            VEL = sv(4:6)
            Cr0 = 1.25_DP               !! 卫星表面反射系数
            SMR0 = 7.0e-3_dp            !! m^2/kg 面质比
            RP0 = 4.5605e-6_dp          !! N/m^2 N = kg*m/s^2 一个AU处的太阳辐射压强度
            SVS = SVDE(:,11)            !! 注意太阳在 SVDE 的编号
            if(present(Cr)) Cr0=Cr
            if(present(SMR)) SMR0 = SMR
            if(present(RP)) RP0 = RP
            SRP = Cr0*SMR0*RP0*(AU**2/norm2(POS-SVS(1:3))**2)*(POS-SVS(1:3))/norm2(POS-SVS(1:3))
            SRP = SRP*1.0e-3_dp         !! 从单位上可见剩余一个 m 没有消去，转为 km 单位
            dsv(1:3) = VEL
            dsv(4:6) = SRP
        end subroutine f_SRP
!!-----------------------------------------------------------------------------------------------
!! 计算后牛顿效应，考虑日地月的一级引力效应
!!
        subroutine f_PNE(SVDE,sv,dsv)
            implicit none
            REAL(DP),intent(in)     ::SVDE(6,11),sv(6)
            real(dp),intent(out)    ::dsv(6)

            integer                 :: i
            real(dp)                :: POS(3),VEL(3),POS0(3),VEL0(3)
            real(dp)                :: F7(3,11)

            POS = sv(1:3)
            VEL = sv(4:6)
        
            do i=1,11
                if(i==3 .or.i==10.or.i==11)then
                    POS0 = pos-SVDE(1:3,i)
                    VEL0 = vel-SVDE(4:6,i)
                    F7(:,i) = GM(i)/c**2/norm2(POS0)**2*((4*GM(i)/norm2(POS0)-norm2(VEL0)**2)&
                            *POS0/norm2(POS0)+4*dot_product(VEL0,POS0/norm2(POS0))*VEL0)
                else
                    F7(:,i)=0.0_dp
                endif
            enddo
            do i=1,3
                dsv(i+3) = sum(F7(i,:))
            enddo
            dsv(1:3) = VEL
        end subroutine f_PNE
!! ===============================================
!! the partial of right function of J2000 equation of motion: 
!! auther:shi zhipeng
!!
!! update history: 2024.6.21 created
!!                 2024.11.18 delete the HFEM part
!! ref : 《航天器动力学引论》
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system[HFEM model will be computed by finite-difference]
!! tMJD       I   MJD时间
!! sv         I   the state vector 
!! DF        O   the differential vector of state vector
!!
    subroutine computeDF_J2000(model,tMJD,sv,DF)
        implicit none

        real(dp),intent(in)::tMJD
        character(len=*),intent(in)::model
        real(dp),dimension(6),intent(in):: sv
        real(dp),dimension(6,6),intent(out)::DF

        !! parameter
        integer::i,j,unitm(3,3)
        real(dp)::SVM(6),SVS(6),dist,DF0(3,3),DF2m(3,3),DF2s(3,3),POS(3),zeros(3,3)
        real(dp)::x,y,z,matrx1(3,3),matrx2(3,3),matrx3(3,3),delta(3),deltax,deltay,deltaz,ldt
        real(dp)::DFF(3,3)! ,SVDE(6,10),matrix(3,3,10),DF2(3,3,10)
        forall(i=1:3,j=1:3)zeros(i,j) = 0.0_dp

        unitm =  reshape((/1,0,0,0,1,0,0,0,1/),shape(unitm))
        POS = sv(1:3)
        x = sv(1)
        y = sv(2)
        z = sv(3)

        call DE440(tMJD,10,3,SVM)
        call DE440(tMJD,11,3,SVS)
        dist = norm2(POS)
        
        matrx1 = reshape((/x**2,x*y,x*z,x*y,y**2,y*z,x*z,y*z,z**2/),shape(matrx1))
        DF0 = -GM(3)/dist**3*(unitm-3/dist**2*matrx1)

        
        DFF = 0.0_dp
        SELECT case(model)
        case('EMS')
        !! MOON
        delta = POS-SVM(1:3)
        deltax = delta(1)
        deltay = delta(2)
        deltaz = delta(3)
        ldt = norm2(delta)
        matrx2 = reshape((/3*deltax**2/ldt**2-1.0_dp,3*deltax*deltay/ldt**2,3*deltax*deltaz/ldt**2,&
                           3*deltax*deltay/ldt**2,3*deltay**2/ldt**2-1.0_dp,3*deltay*deltaz/ldt**2,&
                           3*deltax*deltaz/ldt**2,3*deltay*deltaz/ldt**2,3*deltaz**2/ldt**2-1.0_dp/),shape(matrx2))
        DF2m = GM(11)/ldt**3*matrx2
        !! SUN
        delta = POS-SVS(1:3)
        deltax = delta(1)
        deltay = delta(2)
        deltaz = delta(3)
        ldt = norm2(delta)
        matrx3 = reshape((/3*deltax**2/ldt**2-1.0_dp,3*deltax*deltay/ldt**2,3*deltax*deltaz/ldt**2,&
                           3*deltax*deltay/ldt**2,3*deltay**2/ldt**2-1.0_dp,3*deltay*deltaz/ldt**2,&
                           3*deltax*deltaz/ldt**2,3*deltay*deltaz/ldt**2,3*deltaz**2/ldt**2-1.0_dp/),shape(matrx3))
        DF2s = GM(10)/ldt**3*matrx3
        DFF = DF2m+DF2s
        case('HFEM')
            print*,'wrong here'
        ! j=1
        ! do i = 1,11
        !     if(i/=3)then
        !     call DE440(tMJD,i,3,SVDE(:,j))
        !     delta = POS-SVDE(1:3,j)
        !     deltax = delta(1)
        !     deltay = delta(2)
        !     deltaz = delta(3)
        !     ldt = norm2(delta)
        !     matrix(:,:,j) = reshape((/3*deltax**2/ldt**2-1.0_dp,3*deltax*deltay/ldt**2,3*deltax*deltaz/ldt**2,&
        !                    3*deltax*deltay/ldt**2,3*deltay**2/ldt**2-1.0_dp,3*deltay*deltaz/ldt**2,&
        !                    3*deltax*deltaz/ldt**2,3*deltay*deltaz/ldt**2,3*deltaz**2/ldt**2-1.0_dp/),shape(matrix(:,:,j)))
        !     DF2(:,:,j) = GM(i)/ldt**3*matrix(:,:,j)
        !     j = j+1
        !     endif
        ! enddo
        ! do i=1,10
        !     DFF = DFF+DF2(:,:,i)
        ! enddo
        end SELECT
        

        forall(i=1:3,j=1:3)DF(i,j) = zeros(i,j)
        forall(i=1:3,j=4:6)DF(i,j) = unitm(i,j-3)
        forall(i=4:6,j=1:3)DF(i,j) = DF0(i-3,j)+DFF(I-3,J)
        forall(i=4:6,j=4:6)DF(i,j) = zeros(i-3,j-3)

        end subroutine computeDF_J2000


!!========================================================================
!! the integrator of rkf78
!! auther:shi zhipeng
!!
!! all units is in normalized units
!!
!! update history: 2024.6.24 created
!!
!! Variable  I/O        Description
!! --------  ---        --------------------------------------------------
!!  model     I     the decided model of dynamic system
!!   t        I         the tMJD
!!   x        I         the state vector [allocatable]
!!   h        I         the step length
!!   xf       O         the next step of state vecot
!!  terr      O         the truncation error
!!
!! 
        subroutine rkf78(model,t,x,h,xf,terr)

            implicit none
            character(len=*),intent(in)              :: model
            real(dp),intent(in)                      :: t
            real(dp),intent(in)                      :: x(:)
            real(dp),intent(in)                      :: h
            real(dp),allocatable,intent(out)         :: xf(:)
            real(dp),allocatable,intent(out)                     :: terr(:)

            integer                                  :: dims
        
            real(dp),parameter :: a1  = 2.0_dp/27.0_dp
            real(dp),parameter :: a2  = 1.0_dp/9.0_dp
            real(dp),parameter :: a3  = 1.0_dp/6.0_dp
            real(dp),parameter :: a4  = 5.0_dp/12.0_dp
            real(dp),parameter :: a5  = 1.0_dp/2.0_dp
            real(dp),parameter :: a6  = 5.0_dp/6.0_dp
            real(dp),parameter :: a7  = 1.0_dp/6.0_dp
            real(dp),parameter :: a8  = 2.0_dp/3.0_dp
            real(dp),parameter :: a9  = 1.0_dp/3.0_dp
            !real(dp),parameter :: a10 = 1.0_dp
            !real(dp),parameter :: a12 = 1.0_dp
        
            real(dp),parameter :: b10  = 2.0_dp/27.0_dp
            real(dp),parameter :: b20  = 1.0_dp/36.0_dp
            real(dp),parameter :: b21  = 1.0_dp/12.0_dp
            real(dp),parameter :: b30  = 1.0_dp/24.0_dp
            real(dp),parameter :: b32  = 1.0_dp/8.0_dp
            real(dp),parameter :: b40  = 5.0_dp/12.0_dp
            real(dp),parameter :: b42  = -25.0_dp/16.0_dp
            real(dp),parameter :: b43  = 25.0_dp/16.0_dp
            real(dp),parameter :: b50  = 1.0_dp/20.0_dp
            real(dp),parameter :: b53  = 1.0_dp/4.0_dp
            real(dp),parameter :: b54  = 1.0_dp/5.0_dp
            real(dp),parameter :: b60  = -25.0_dp/108.0_dp
            real(dp),parameter :: b63  = 125.0_dp/108.0_dp
            real(dp),parameter :: b64  = -65.0_dp/27.0_dp
            real(dp),parameter :: b65  = 125.0_dp/54.0_dp
            real(dp),parameter :: b70  = 31.0_dp/300.0_dp
            real(dp),parameter :: b74  = 61.0_dp/225.0_dp
            real(dp),parameter :: b75  = -2.0_dp/9.0_dp
            real(dp),parameter :: b76  = 13.0_dp/900.0_dp
            real(dp),parameter :: b80  = 2.0_dp
            real(dp),parameter :: b83  = -53.0_dp/6.0_dp
            real(dp),parameter :: b84  = 704.0_dp/45.0_dp
            real(dp),parameter :: b85  = -107.0_dp/9.0_dp
            real(dp),parameter :: b86  = 67.0_dp/90.0_dp
            real(dp),parameter :: b87  = 3.0_dp
            real(dp),parameter :: b90  = -91.0_dp/108.0_dp
            real(dp),parameter :: b93  = 23.0_dp/108.0_dp
            real(dp),parameter :: b94  = -976.0_dp/135.0_dp
            real(dp),parameter :: b95  = 311.0_dp/54.0_dp
            real(dp),parameter :: b96  = -19.0_dp/60.0_dp
            real(dp),parameter :: b97  = 17.0_dp/6.0_dp
            real(dp),parameter :: b98  = -1.0_dp/12.0_dp
            real(dp),parameter :: b100 = 2383.0_dp/4100.0_dp
            real(dp),parameter :: b103 = -341.0_dp/164.0_dp
            real(dp),parameter :: b104 = 4496.0_dp/1025.0_dp
            real(dp),parameter :: b105 = -301.0_dp/82.0_dp
            real(dp),parameter :: b106 = 2133.0_dp/4100.0_dp
            real(dp),parameter :: b107 = 45.0_dp/82.0_dp
            real(dp),parameter :: b108 = 45.0_dp/164.0_dp
            real(dp),parameter :: b109 = 18.0_dp/41.0_dp
            real(dp),parameter :: b110 = 3.0_dp/205.0_dp
            real(dp),parameter :: b115 = -6.0_dp/41.0_dp
            real(dp),parameter :: b116 = -3.0_dp/205.0_dp
            real(dp),parameter :: b117 = -3.0_dp/41.0_dp
            real(dp),parameter :: b118 = 3.0_dp/41.0_dp
            real(dp),parameter :: b119 = 6.0_dp/41.0_dp
            real(dp),parameter :: b120 = -1777.0_dp/4100.0_dp
            real(dp),parameter :: b123 = -341.0_dp/164.0_dp
            real(dp),parameter :: b124 = 4496.0_dp/1025.0_dp
            real(dp),parameter :: b125 = -289.0_dp/82.0_dp
            real(dp),parameter :: b126 = 2193.0_dp/4100.0_dp
            real(dp),parameter :: b127 = 51.0_dp/82.0_dp
            real(dp),parameter :: b128 = 33.0_dp/164.0_dp
            real(dp),parameter :: b129 = 12.0_dp/41.0_dp
            !real(dp),parameter :: b1211 = 1.0_dp
        
            real(dp),parameter :: c5  = 34.0_dp/105.0_dp
            real(dp),parameter :: c6  = 9.0_dp/35.0_dp
            real(dp),parameter :: c7  = 9.0_dp/35.0_dp
            real(dp),parameter :: c8  = 9.0_dp/280.0_dp
            real(dp),parameter :: c9  = 9.0_dp/280.0_dp
            real(dp),parameter :: c11 = 41.0_dp/840.0_dp
            real(dp),parameter :: c12 = 41.0_dp/840.0_dp
        
            real(dp),allocatable :: f0(:),f1(:),f2(:),f3(:),f4(:),f5(:),&
                                f6(:),f7(:),f8(:),f9(:),f10(:),f11(:),f12(:)
            real(dp)             :: hs   !! second unit of h

            if (h==0) then
                xf = x
                terr = 0.0_dp
                return
            end if

            hs = h*86400.0_dp
            dims = size(x,dim=1)

            allocate(xf(dims))
            allocate(f0(dims))
            allocate(f1(dims))
            allocate(f2(dims))
            allocate(f3(dims))
            allocate(f4(dims))
            allocate(f5(dims))
            allocate(f6(dims))
            allocate(f7(dims))
            allocate(f8(dims))
            allocate(f9(dims))
            allocate(f10(dims))
            allocate(f11(dims))
            allocate(f12(dims))
            allocate(terr(dims))
            
            !! in situation of crtbp the f dont have t
            
            call f_EPH(model,t,x,f0)
            call f_EPH(model,t+h*a1,x+f0*b10*hs,f1)
            call f_EPH(model,t+h*a2,x+(f0*b20+f1*b21)*hs,f2)
            call f_EPH(model,t+h*a3,x+(f0*b30+f2*b32)*hs,f3)
            call f_EPH(model,t+h*a4,x+(f0*b40+f2*b42+f3*b43)*hs,f4)
            call f_EPH(model,t+h*a5,x+(f0*b50+f3*b53+f4*b54)*hs,f5)
            call f_EPH(model,t+h*a6,x+(f0*b60+f3*b63+f4*b64+f5*b65)*hs,f6)
            call f_EPH(model,t+h*a7,x+(f0*b70+f4*b74+f5*b75+f6*b76)*hs,f7)
            call f_EPH(model,t+h*a8,x+(f0*b80+f3*b83+f4*b84+f5*b85+f6*b86+&
                        f7*b87)*hs,f8)
            call f_EPH(model,t+h*a9,x+(f0*b90+f3*b93+f4*b94+f5*b95+f6*b96+&
                        f7*b97+f8*b98)*hs,f9)
            call f_EPH(model,t+h,x+(f0*b100+f3*b103+f4*b104+f5*b105+&
                        f6*b106+f7*b107+f8*b108+f9*b109)*hs,f10)
            call f_EPH(model,t,x+(f0*b110+f5*b115+f6*b116+f7*b117+f8*b118+&
                        f9*b119)*hs,f11)
            call f_EPH(model,t+h,x+(f0*b120+f3*b123+f4*b124+f5*b125+f6*b126+&
                        f7*b127+f8*b128+f9*b129+f11)*hs,f12)


            xf = x + hs*(f5*c5+f6*c6+f7*c7+f8*c8+f9*c9+f11*c11+f12*c12)
            ! write(*,*) "xf",xf
            ! print*,size(f0,1)
            terr = (41.0_dp/840.0_dp)*(f0+f10-f11-f12)*abs(hs)
            ! write(*,*) "terr",terr
            end subroutine rkf78

!!===========================================================================
!! propagateOrbit_crtbp part
!! auther:shi zhipeng
!!
!! all units is in normalized units
!!
!! update history: 2024.4.17 created
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system
!!   t        I      tMJD
!!   svn      I      the initial state vector [6 or 42]
!!  svnd      O      the differential state vector [6 or 42]

            subroutine f_EPH(model,t,svn,svnd)
                implicit none
        
                real(dp),intent(in)                 :: t
                character(len=*),intent(IN)         :: model
                real(dp),intent(in)                 :: svn(:)
                real(dp),allocatable,intent(out)    :: svnd(:)
        
                real(dp),dimension(6)               :: sv          !! the postion part of sv
                real(dp),dimension(6)               :: dsv          
                real(dp),dimension(6,6)             :: DF,dSTM,STM           
        
                integer                             :: dims
        
                dims = size(svn,1)
                allocate(svnd(dims))
        
                if (size(svn,1)==6) then
                    call computeF_J2000(model,t,svn,dsv)
                    svnd = dsv
                else if (size(svn,1)==42 .and. model /='HFEM') then
                    sv = svn(1:6)
                    STM = reshape(svn(7:42),[6,6])
                    call computeF_J2000(model,t,sv,dsv)
                    call computeDF_J2000(model,t,sv,DF)
                    dSTM = matmul(DF,STM)
        
                    svnd(1:6) = dsv
                    svnd(7:42) = reshape(dSTM,[36])
                else
                    print *, 'the dimension of sv is wrong'
                end if
        
                end subroutine f_EPH
     
!!===========================================================================
!! propagateOrbit_EPH part of MIDDLE transitions for HFEM and EMS 's different STM computing method
!! auther:shi zhipeng
!!
!!
!! update history: 2024.4.17 created
!!                 2024.6.17 created sxh to Avoid continuous changes in h
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system
!!   t        I      tMJD
!!   sv0      I      the initial state vector [6]
!!   tf       I      the propagate time tf[days]
!! sv_nd      O     the state vector in every points of integration[6xN]
!!  tol       I     the truncation error
!!   t        O     the time of every integration points[N]
!!  STM     O[op]   the state transition matrix [optional] [6x6xN][EMS]
!!
    subroutine propagateOrbit_EPH_MIDDLE(model,tMJD,sv0,tf,tol,sv_nd,t,STM)
        implicit none

        character(len=*),intent(in)                         :: model
        real(dp),intent(in)                                 :: tMJD
        real(dp),intent(in)                                 :: sv0(6)
        real(dp),intent(in)                                 :: tf
        real(dp),intent(in)                                 :: tol

        real(dp),allocatable,intent(out)                    :: sv_nd(:,:)
        real(dp),optional,allocatable,intent(out)           :: STM(:,:,:)
        real(dp),allocatable,intent(out)                    :: t(:)

        real(dp),allocatable                                :: svn(:),svnn(:,:),terr(:)        !! the vector used in dieration for transition 
        integer                                             :: I (6,6)              !! the I matrix
        integer                                             :: k,j,sxh
        real(dp)                                            :: h!,terr(6)   !! the initial step length and for Fine processing and tol in precedure
        real(dp),parameter                                  :: hmin = 1.0e-6_dp,hmax = 12.0/24.0_dp !! the max and min stepszie:1h and 0.0864s

        !! derive a 1 matrix to the initial matrix of stm
        forall(k=1:6,j=1:6,j==k)  I(k,j) = 1
        forall(k=1:6,j=1:6,j/=k)  I(k,j) = 0
        sxh =0
            h=tf/abs(tf)*20.0_dp/86400_dp  !! 初始20s的步长
            allocate(sv_nd(6,int(abs(tf)/hmin)+2))
            allocate(t(int(abs(tf)/hmin)+2))
        
            if(present(STM))then
                allocate(svn(42))
                allocate(svnn(42,int(abs(tf)/hmin)+2))
                allocate(STM(6,6,int(abs(tf)/hmin)+2))
                
                svn(1:6) = sv0
                svn(7:42) = reshape(I,[36])

                sv_nd(:,1) = sv0
                STM(:,:,1) = I
               
                svnn(:,1) = svn
                t(1) =0.0_dp
                k=2
                do while(k<abs(tf)/hmin)
                    call rkf78(model,tMJD+t(k-1),svnn(:,k-1),h,svn,terr)
                    sv_nd(:,k) = svn(1:6)
                    STM(:,:,k) = reshape(svn(7:42),[6,6])
                    svnn(:,k) = svn 
                    t(k) = t(k-1)+h
                    if(abs(t(k)-tf)<1e-20_dp)then
                        exit
                    endif
                    ! print*,svn
                    if(norm2(terr(1:6))>tol)then
                        ! if(norm2(terr)>tol)then
                        h=h/2
                        k=k-1
                        sxh=sxh+1
                    elseif(norm2(terr(1:6))<tol*1.0e-15_dp)then
                        ! elseif(norm2(terr)<tol*1.0e-15_dp)then
                        h=h*2
                        k=k-1
                        sxh=sxh+1
                    else
                        sxh =0
                    endif
                    if(sxh>3)then
                        h = h/3
                        sxh=0
                    endif
                    if(tf>=0 .and.t(k)+h>tf .and. t(k)<tf)then
                        h=tf-t(k)
                    elseif(tf>=0 .and.t(k)+h>tf .and. t(k)>tf)then  !! 解决固定步长输出时初始步长大于tf时死循环
                        h = tf/2
                        k =k-1
                    elseif(tf<0 .and. t(k)+h<tf .and. t(k)>tf)then
                        h =tf-t(k)
                    elseif(tf<0 .and. t(k)+h<tf .and. t(k)<tf)then  !! 解决固定步长输出时初始步长大于tf时死循环
                        h = tf/2
                        k =k-1
                    endif
                    ! write(*,*)'步长(s),步数,时间,误差',h*86400,k,t(k),norm2(terr)
                    k=k+1
                    !write(*,*) sxh
                enddo
                sv_nd = sv_nd(:,1:k)
                STM = STM(:,:,1:k)
                t = t(1:k)+tMJD
                deallocate(svnn)
            else
                allocate(svn(6))
                svn(1:6) = sv0
                sv_nd(:,1) = sv0
                t(1) =0.0_dp
                k=2
                do while(k<abs(tf)/hmin)
                    call rkf78(model,tMJD+t(k-1),sv_nd(:,k-1),h,svn,terr)
                    sv_nd(:,k) = svn(1:6)
                    t(k) = t(k-1)+h
                    if(abs(t(k)-tf)<1e-20_dp)then
                        exit
                    endif
                    if(norm2(terr)>tol)then! (1:6)
                        h=h/2
                        k=k-1
                        sxh=sxh+1
                    elseif(norm2(terr)<tol*1.0e-15_dp)then!(1:6)
                        h=h*2
                        k=k-1
                        sxh=sxh+1
                    else
                        sxh =0
                    endif
                    if(sxh>3)then
                        h = h/3
                        sxh=0
                    endif
                    if(tf>=0 .and.t(k)+h>tf .and. t(k)<tf)then
                        h=tf-t(k)
                    elseif(tf>=0 .and.t(k)+h>tf .and. t(k)>tf)then  !! 解决固定步长输出时初始步长大于tf时死循环
                        h = tf/2
                        k =k-1
                    elseif(tf<0 .and. t(k)+h<tf .and. t(k)>tf)then
                        h =tf-t(k)
                    elseif(tf<0 .and. t(k)+h<tf .and. t(k)<tf)then  !! 解决固定步长输出时初始步长大于tf时死循环
                        h = tf/2
                        k =k-1
                    endif
                    ! write(*,*)'h=',h,k,t(k),norm2(terr),sv_nd(:,k)
                    ! open(unit=1,file="C:\Users\28793\Desktop\h7.txt")
                    ! write(1,*) h,k,t(k),norm2(terr),sv_nd(:,k)
                    
                    k=k+1
                enddo
            sv_nd = sv_nd(:,1:k)
            t = t(1:k)+tMJD
            endif
            deallocate(svn)
            ! close(1)
        end subroutine propagateOrbit_EPH_MIDDLE  

!!=============================================================================
!! use the finite difference method to compute STM
!! auther : shi zhipeng
!!
!! created: 2024.11.18
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system
!!   t        I      tMJD
!!   sv0      I      the initial state vector [6]
!!   tf       I      the propagate time tf[days]
!!  tol       I     the truncation error
!! sv_nd      O     the state vector in every points of integration[6xN]
!!   t        O     the time of every integration points[N]
!!  STM       O   the state transition matrix  [6x6]
!! 
        subroutine finite_diff(model,tMJD,sv0,tf,tol,SV,t,STM)
            implicit none
            character(len=*),intent(in)                         :: model
            real(dp),intent(in)                                 :: tMJD
            real(dp),intent(in)                                 :: sv0(6)
            real(dp),intent(in)                                 :: tf
            real(dp),intent(in)                                 :: tol
    
            real(dp),allocatable,intent(out)                    :: SV(:,:)
            real(dp),intent(out)                                :: STM(6,6)
            real(dp),allocatable,intent(out)                    :: t(:)

            integer                                             :: I(6,6),j,k
            real(dp)                                            :: deltar,deltav
            real(dp),allocatable                                :: SVV(:,:),tt(:)
            I = 0
            forall(j=1:6,k=1:6,j==k)I(j,k) =1
            deltar = 1.0e-5_dp      !! 1cm
            deltav = 1.0e-5_dp      !! 10^-5 km/s    
            call propagateOrbit_EPH_middle(model,tMJD,sv0,tf,tol,SV,t)

            do j = 1,3
                call propagateOrbit_EPH_middle(model,tMJD,sv0+deltar*I(:,j),tf,tol,SVV,tt)
                STM(:,j) = (SVV(:,size(SVV,2))-SV(:,size(SV,2)))/deltar
            enddo
            deallocate(SVV,tt)
            do j =4,6
                call propagateOrbit_EPH_middle(model,tMJD,sv0+deltav*I(:,j),tf,tol,SVV,tt)
                STM(:,j) = (SVV(:,size(SVV,2))-SV(:,size(SV,2)))/deltar
            enddo
            deallocate(SVV,tt)
        end subroutine finite_diff
!!===========================================================================
!! propagateOrbit_EPH part 
!! auther:shi zhipeng
!!
!!
!! update history: 2024.11.18 created
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system
!!   t        I      tMJD
!!   sv0      I      the initial state vector [6]
!!   tf       I      the propagate time tf[days]
!!    SV      O     the state vector in every points of integration[6xN]
!!  tol       I     the truncation error
!!   t        O     the time of every integration points[N]
!!  STM     O[op]   the state transition matrix [optional] [6x6]
!! 
        subroutine propagateOrbit_EPH(model,tMJD,sv0,tf,tol,SV,t,STM)
            implicit none

            character(len=*),intent(in)                         :: model
            real(dp),intent(in)                                 :: tMJD
            real(dp),intent(in)                                 :: sv0(6)
            real(dp),intent(in)                                 :: tf
            real(dp),intent(in)                                 :: tol
    
            real(dp),allocatable,intent(out)                    :: SV(:,:)
            real(dp),optional,intent(out)                       :: STM(6,6)
            real(dp),allocatable,intent(out)                    :: t(:)

            real(dp),allocatable                                :: STMM(:,:,:)


            SELECT case(model)
            case('EMS')
                if(present(STM))then
                ! call propagateOrbit_EPH_middle(model,tMJD,sv0,tf,tol,SV,t,STMM)
                ! STM = STMM(:,:,size(STMM,3))
                call finite_diff(model,tMJD,sv0,tf,tol,SV,t,STM)
                else
                call propagateOrbit_EPH_middle(model,tMJD,sv0,tf,tol,SV,t)
                endif
            case('HFEM')
                IF(present(STM))then
                call finite_diff(model,tMJD,sv0,tf,tol,SV,t,STM)
                else
                call propagateOrbit_EPH_middle(model,tMJD,sv0,tf,tol,SV,t)
                endif
            end SELECT
            if(allocated(STMM)) deallocate(STMM)
            
        end subroutine propagateOrbit_EPH
!!!!===========================================================================
!! propagateOrbit_J2000 part in a fixed step length to output a Orbit
!! auther:shi zhipeng
!!
!! the long time propagate may lead to more inaccuracies
!!
!! update history: 2024.9.20 created
!!
!! Variable  I/O  Description
!! --------  ---  --------------------------------------------------
!!  model     I     the decided model of dynamic system
!!   tMJD     I     begin tMJD
!!   sv0      I      the initial state vector [6]
!!   tf       I      the propagate time tf[days]
!!   step     I     the step length[mins]
!! sv_nd      O     the state vector in every points of integration[6xN]
!!  tol       I     the truncation error
!!   t        O     the time of every integration points[N]
!!
        subroutine propagateOrbit_EPH_fixedstep(model,tMJD,sv0,tf,steplength,tol,sv_nd,t)
            implicit none

            character(len=*),intent(in)                         :: model
            real(dp),intent(in)                                 :: sv0(6)
            real(dp),intent(in)                                 :: tf,tMJD
            real(dp),intent(in)                                 :: tol
            real(dp),intent(in)                                 :: steplength
    
            real(dp),allocatable,intent(out)                    :: sv_nd(:,:)
            real(dp),allocatable,intent(out)                    :: t(:)

            real(dp),allocatable                                :: svn(:,:),tn(:)
            integer                                             :: i,k
            real(dp)                                            :: hn,step

            ! real(dp)                                            :: unitim,ET,SVM(6)

            step = steplength/24.0_dp/60.0_dp
            k = int(tf/step)
            hn = mod(tf,step)
            ! print*,k,hn
            if (hn<1.0e-10_dp)then
                allocate(sv_nd(6,k+1),t(k+1))
            else
                allocate(sv_nd(6,k+2),t(k+2))
            endif
            sv_nd(:,1) = sv0
            t(1) = 0.0_dp+tMJD
            do i=1,k
                ! print*,i,k
                call propagateOrbit_EPH(model,t(i),sv_nd(:,i),step,tol,svn,tn)
                sv_nd(:,i+1) = svn(:,size(svn,2))
                t(i+1) = i*step+tMJD

            enddo

            if(hn>1.0e-10_dp)then
                call propagateOrbit_EPH(model,t(k+1),sv_nd(:,k+1),hn,tol,svn,tn)
                sv_nd(:,k+2) = svn(:,size(svn,2))
                t(k+2) = k*step+tMJD+hn
            endif
            if(allocated(svn)) deallocate(svn)
            if(allocated(tn)) deallocate(tn)
        end subroutine propagateOrbit_EPH_fixedstep


!!---------------------------------------------------------------------------

end module propagateOrbit_EPH_module
!!===========================================================    