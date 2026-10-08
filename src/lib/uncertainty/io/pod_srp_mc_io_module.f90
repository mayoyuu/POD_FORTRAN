!> Strict CSV I/O for fixed SRP Monte Carlo samples and point staging files.
module pod_srp_mc_io_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use pod_global, only: DP
    implicit none
    private

    character(len=*),parameter,public :: SRP_MC_POINTS_SCHEMA='srp_mc_points_v1'
    character(len=*),parameter,public :: SRP_MC_INITIAL_POINTS_SCHEMA= &
        'srp_mc_initial_points_v1'
    character(len=*),parameter :: SAMPLE_HEADER= &
        'sample_id,z1,z2,z3,z4,z5,z6,z7,u1,u2,u3,u4,u5,u6,u7,eta_treatment'
    real(DP),parameter :: SIGMA_LIMIT=3.0_DP,DEFAULT_ETA_SIGMA=0.3_DP

    public :: read_srp_mc_samples,write_mc_point_header,write_mc_point_rows
    public :: write_initial_only_point_header,write_initial_only_point_rows

contains

    subroutine read_srp_mc_samples(path,expected_count,z,u,status,message,expected_eta_sigma)
        character(len=*),intent(in) :: path
        integer,intent(in) :: expected_count
        real(DP),allocatable,intent(out) :: z(:,:),u(:,:)
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        real(DP),intent(in),optional :: expected_eta_sigma

        character(len=4096) :: line
        real(DP) :: eta,sample_eta_sigma
        integer :: input_unit,ios,row,sample_id

        status=0
        message=''
        sample_eta_sigma=DEFAULT_ETA_SIGMA
        if(present(expected_eta_sigma)) sample_eta_sigma=expected_eta_sigma
        if(allocated(z)) deallocate(z)
        if(allocated(u)) deallocate(u)
        if(.not.ieee_is_finite(sample_eta_sigma).or.sample_eta_sigma<0.0_DP) then
            call fail(-13,'expected eta sigma must be finite and nonnegative')
            return
        end if
        if(expected_count<=0) then
            call fail(-1,'expected sample count must be positive')
            return
        end if
        open(newunit=input_unit,file=trim(path),status='old',action='read',iostat=ios)
        if(ios/=0) then
            call fail(-2,'cannot open sample CSV')
            return
        end if
        read(input_unit,'(A)',iostat=ios) line
        if(ios/=0 .or. trim(line)/=SAMPLE_HEADER) then
            close(input_unit)
            call fail(-3,'sample CSV header does not match srp_mc_samples_v1')
            return
        end if

        allocate(z(7,expected_count),u(7,expected_count))
        z=0.0_DP
        u=0.0_DP
        do row=1,expected_count
            read(input_unit,'(A)',iostat=ios) line
            if(ios/=0) then
                close(input_unit)
                call fail(-4,'sample CSV ended before expected row count')
                return
            end if
            if(count_character(trim(line),',')/=15) then
                close(input_unit)
                call fail(-5,'sample row has the wrong number of columns')
                return
            end if
            read(line,*,iostat=ios) sample_id,z(:,row),u(:,row),eta
            if(ios/=0) then
                close(input_unit)
                call fail(-6,'sample row is malformed')
                return
            end if
            if(sample_id/=row) then
                close(input_unit)
                call fail(-7,'sample IDs must be unique and sequential from one')
                return
            end if
            if(any(.not.ieee_is_finite(z(:,row))) .or. &
               any(.not.ieee_is_finite(u(:,row))) .or. .not.ieee_is_finite(eta)) then
                close(input_unit)
                call fail(-8,'sample row contains a nonfinite value')
                return
            end if
            if(any(abs(z(:,row))>SIGMA_LIMIT)) then
                close(input_unit)
                call fail(-9,'sample coordinate exceeds componentwise three sigma limit')
                return
            end if
            if(maxval(abs(u(:,row)-z(:,row)/SIGMA_LIMIT))>1.0e-12_DP) then
                close(input_unit)
                call fail(-10,'sample normalized coordinate is inconsistent with u=z/3')
                return
            end if
            if(abs(eta-sample_eta_sigma*z(7,row))>1.0e-12_DP* &
               max(1.0_DP,abs(eta),abs(sample_eta_sigma*z(7,row)))) then
                close(input_unit)
                call fail(-11,'sample eta is inconsistent with eta=sigma_eta*z7')
                return
            end if
        end do
        read(input_unit,'(A)',iostat=ios) line
        close(input_unit)
        if(ios==0) then
            call fail(-12,'sample CSV has more rows than expected')
            return
        end if

    contains

        subroutine fail(code,text)
            integer,intent(in) :: code
            character(len=*),intent(in) :: text
            status=code
            message=text
            if(allocated(z)) deallocate(z)
            if(allocated(u)) deallocate(u)
        end subroutine fail

    end subroutine read_srp_mc_samples

    subroutine write_mc_point_header(unit_value,status,message)
        integer,intent(in) :: unit_value
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        integer :: ios

        status=0
        message=''
        write(unit_value,'(A)',iostat=ios) 'schema='//SRP_MC_POINTS_SCHEMA
        if(ios==0) write(unit_value,'(A)',iostat=ios) &
            'hour,keyframe_kind,sample_id,'// &
            'treatment_dx_km,treatment_dy_km,treatment_dz_km,'// &
            'treatment_dvx_kms,treatment_dvy_kms,treatment_dvz_kms,'// &
            'control_dx_km,control_dy_km,control_dz_km,'// &
            'control_dvx_kms,control_dvy_kms,control_dvz_kms,'// &
            'delta_dx_km,delta_dy_km,delta_dz_km,'// &
            'delta_dvx_kms,delta_dvy_kms,delta_dvz_kms'
        if(ios/=0) then
            status=-1
            message='failed to write MC point header'
        end if
    end subroutine write_mc_point_header

    subroutine write_mc_point_rows(unit_value,hour,keyframe_kind,sample_ids, &
                                   treatment_errors,control_errors,status,message)
        integer,intent(in) :: unit_value,sample_ids(:)
        real(DP),intent(in) :: hour,treatment_errors(:,:),control_errors(:,:)
        character(len=*),intent(in) :: keyframe_kind
        integer,intent(out) :: status
        character(len=*),intent(out) :: message

        real(DP) :: delta(6)
        integer :: row_count,row,ios

        status=0
        message=''
        row_count=size(sample_ids)
        if(.not.ieee_is_finite(hour) .or. len_trim(keyframe_kind)==0 .or. &
           index(keyframe_kind,',')>0 .or. size(treatment_errors,1)/=6 .or. &
           size(control_errors,1)/=6 .or. size(treatment_errors,2)/=row_count .or. &
           size(control_errors,2)/=row_count .or. row_count<1) then
            status=-1
            message='invalid MC point row dimensions or keyframe metadata'
            return
        end if
        if(any(.not.ieee_is_finite(treatment_errors)) .or. &
           any(.not.ieee_is_finite(control_errors))) then
            status=-2
            message='MC point rows contain nonfinite values'
            return
        end if
        do row=1,row_count
            delta=treatment_errors(:,row)-control_errors(:,row)
            write(unit_value,'(ES25.16E3,",",A,",",I0,18(",",ES25.16E3))', &
                iostat=ios) hour,trim(keyframe_kind),sample_ids(row), &
                treatment_errors(:,row),control_errors(:,row),delta
            if(ios/=0) then
                status=-3
                message='failed to write MC point row'
                return
            end if
        end do
    end subroutine write_mc_point_rows

    subroutine write_initial_only_point_header(unit_value,status,message)
        integer,intent(in) :: unit_value
        integer,intent(out) :: status
        character(len=*),intent(out) :: message
        integer :: ios

        status=0
        message=''
        write(unit_value,'(A)',iostat=ios) 'schema='//SRP_MC_INITIAL_POINTS_SCHEMA
        if(ios==0) write(unit_value,'(A)',iostat=ios) &
            'hour,keyframe_kind,sample_id,'// &
            'x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,'// &
            'dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms'
        if(ios/=0) then
            status=-1
            message='failed to write initial-only MC point header'
        end if
    end subroutine write_initial_only_point_header

    subroutine write_initial_only_point_rows(unit_value,hour,keyframe_kind,sample_ids, &
                                             states,center,status,message)
        integer,intent(in) :: unit_value,sample_ids(:)
        real(DP),intent(in) :: hour,states(:,:),center(:)
        character(len=*),intent(in) :: keyframe_kind
        integer,intent(out) :: status
        character(len=*),intent(out) :: message

        real(DP) :: error_state(6)
        integer :: row_count,row,ios

        status=0
        message=''
        row_count=size(sample_ids)
        if(.not.ieee_is_finite(hour).or.len_trim(keyframe_kind)==0.or. &
           index(keyframe_kind,',')>0.or.size(states,1)/=6.or.size(center)/=6.or. &
           size(states,2)/=row_count.or.row_count<1) then
            status=-1
            message='invalid initial-only MC point dimensions or keyframe metadata'
            return
        end if
        if(any(.not.ieee_is_finite(states)).or.any(.not.ieee_is_finite(center))) then
            status=-2
            message='initial-only MC point rows contain nonfinite values'
            return
        end if
        do row=1,row_count
            error_state=states(:,row)-center
            write(unit_value,'(ES25.16E3,",",A,",",I0,12(",",ES25.16E3))', &
                iostat=ios) hour,trim(keyframe_kind),sample_ids(row), &
                states(:,row),error_state
            if(ios/=0) then
                status=-3
                message='failed to write initial-only MC point row'
                return
            end if
        end do
    end subroutine write_initial_only_point_rows
    pure integer function count_character(text,needle) result(count_value)
        character(len=*),intent(in) :: text
        character,intent(in) :: needle
        integer :: index_value

        count_value=0
        do index_value=1,len_trim(text)
            if(text(index_value:index_value)==needle) count_value=count_value+1
        end do
    end function count_character

end module pod_srp_mc_io_module
