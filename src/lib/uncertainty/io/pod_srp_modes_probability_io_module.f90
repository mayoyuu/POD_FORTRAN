!> Opt-in probability IO. Population statistics never use the visualization subset.
module pod_srp_modes_probability_io_module
    use pod_global, only: DP
    use pod_srp_mc_statistics_module, only: srp_mc_group_stats_type,compute_group_stats
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    character(len=*),parameter,public :: PROBABILITY_SAMPLE_HEADER= &
        'sample_id,z1,z2,z3,z4,z5,z6,q1,q2,q3,q4,q5,q6,q7,q8,s,u,xi_E,xi_M'
    character(len=*),parameter :: POINT_HEADER= &
        'hour,epoch_et,keyframe_kind,sample_id,patch_id,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms,'// &
        'dx_km,dy_km,dz_km,dvx_kms,dvy_kms,dvz_kms'
    type,public :: probability_population_type
        real(DP),allocatable :: q(:,:),modes(:,:)
    end type
    type,public :: probability_output_type
        integer :: stats_unit=0,visual_unit=0,key_unit=0
        integer,allocatable :: visual_ids(:),all_ids(:),prefix_counts(:),previous_owners(:)
        real(DP),allocatable :: previous_states(:,:)
        real(DP) :: previous_hour=-1,previous_epoch=0,previous_reference(6)=0
        real(DP) :: first_failure_hour=-1,first_unconfirmed_hour=-1
        character(len=256) :: previous_kind=''
        logical :: seen_split=.false.,seen_depth=.false.,seen_failure=.false.
    end type
    public :: load_probability_samples,probability_choose_ids,probability_open,probability_write_hour,probability_close
contains
    subroutine probability_choose_ids(n,requested,ids,status)
        integer,intent(in) :: n,requested
        integer,allocatable,intent(out) :: ids(:)
        integer,intent(out) :: status
        integer :: j,v
        status=-1
        if(min(n,requested)<1) return
        v=min(n,requested); allocate(ids(v)); ids=1
        if(v>1) then
            do j=1,v
                ids(j)=1+nint(real(j-1,DP)*real(n-1,DP)/real(v-1,DP))
            end do
        end if
        status=0
    end subroutine

    subroutine load_probability_samples(path,n,population,status)
        character(len=*),intent(in) :: path
        integer,intent(in) :: n
        type(probability_population_type),intent(out) :: population
        integer,intent(out) :: status
        character(len=8192) :: line
        character(len=100) :: token
        real(DP) :: row(19)
        integer :: unit,io,i,j,start,finish,k
        status=-1
        if(n<2) return
        open(newunit=unit,file=path,status='old',action='read',iostat=io)
        if(io/=0) return
        read(unit,'(A)',iostat=io) line
        if(io==0) then
            if(trim(line)/='schema=srp_modes_probability_samples_v1') io=1
        end if
        if(io==0) read(unit,'(A)',iostat=io) line
        if(io==0) then
            if(trim(line)/=PROBABILITY_SAMPLE_HEADER) io=1
        end if
        if(io/=0) then
            close(unit); return
        end if
        allocate(population%q(8,n),population%modes(2,n))
        do i=1,n
            read(unit,'(A)',iostat=io) line
            if(io/=0) exit
            if(count([(line(k:k)==',',k=1,len_trim(line))])/=18) then
                io=1; exit
            end if
            start=1
            do j=1,19
                finish=index(line(start:),',')+start-2
                if(j==19) finish=len_trim(line)
                if(finish<start.or.finish-start+1>len(token)) then
                    io=1; exit
                end if
                token=trim(adjustl(line(start:finish)))
                if(len_trim(token)==0.or.verify(trim(token),'0123456789.+-EeDd')/=0) then
                    io=1; exit
                end if
                read(token,*,iostat=io) row(j)
                if(io/=0) exit
                start=finish+2
            end do
            if(io/=0) exit
            if(.not.all(ieee_is_finite(row))) then
                io=1; exit
            end if
            if(row(1)/=real(i,DP).or.any(abs(row(2:7))>3).or.any(abs(row(8:15))>1).or. &
                any(row(16:17)<0).or.any(row(16:17)>1)) then
                io=1; exit
            end if
            if(any(abs(row(8:13)-row(2:7)/3)>2e-14_DP).or. &
                any(abs(row(14:15)-(2*row(16:17)-1))>2e-14_DP).or. &
                abs(row(18)-row(16))>2e-14_DP.or.abs(row(19)-(1-row(16))*row(17))>2e-14_DP) then
                io=1; exit
            end if
            population%q(:,i)=row(8:15); population%modes(:,i)=row(16:17)
        end do
        if(io==0) then
            read(unit,'(A)',iostat=io) line
            if(io<0) status=0
        end if
        close(unit)
        if(status/=0) deallocate(population%q,population%modes)
    end subroutine

    subroutine probability_open(output,prefix,population,requested,status)
        type(probability_output_type),intent(out) :: output
        character(len=*),intent(in) :: prefix
        type(probability_population_type),intent(in) :: population
        integer,intent(in) :: requested
        integer,intent(out) :: status
        integer :: n,i,j,k,standard(3),unit
        character(len=:),allocatable :: header
        n=size(population%q,2)
        call probability_choose_ids(n,requested,output%visual_ids,status)
        if(status/=0.or.n<2) then
            status=-1; return
        end if
        allocate(output%all_ids(n),output%previous_owners(n),output%previous_states(6,n))
        output%all_ids=[(i,i=1,n)]
        standard=[10000,20000,50000]
        output%prefix_counts=[pack(standard,standard<n),n]
        open(newunit=output%stats_unit,file=prefix//'_probability_hourly_stats.csv',status='new',iostat=status)
        if(status/=0) return
        open(newunit=output%visual_unit,file=prefix//'_probability_visualization.csv',status='new',iostat=status)
        if(status/=0) return
        open(newunit=output%key_unit,file=prefix//'_probability_keyframes.csv',status='new',iostat=status)
        if(status/=0) return
        do i=1,2
            unit=output%visual_unit
            if(i==2) unit=output%key_unit
            write(unit,'(A)') 'schema=srp_modes_probability_points_v1'
            write(unit,'(A)') POINT_HEADER
        end do
        header='hour,epoch_et,statistics_sample_count,saved_visualization_count'
        call header_vector(header,'reference',6); call header_vector(header,'mean_state',6)
        call header_vector(header,'mean_error',6)
        do i=1,6
            do j=1,6
                header=header//',cov_'//digit(i)//digit(j)
            end do
        end do
        call header_vector(header,'std',6); call header_vector(header,'min',6); call header_vector(header,'max',6)
        call header_vector(header,'skewness',6); call header_vector(header,'excess_kurtosis',6)
        call header_vector(header,'skewness_valid',6); call header_vector(header,'kurtosis_valid',6)
        call header_vector(header,'position_axis_rms',3); call header_vector(header,'velocity_axis_rms',3)
        header=header//',position_norm_rms,velocity_norm_rms'
        do k=1,2
            if(k==1) then
                header=header//',position_p50,position_p95,position_p99,position_p999,position_max,position_worst_id'
            else
                header=header//',velocity_p50,velocity_p95,velocity_p99,velocity_p999,velocity_max,velocity_worst_id'
            end if
        end do
        header=header//',patches,depth_limited,original_probe_pos_km,original_probe_vel_kms,'// &
            'probability_probe_pos_km,probability_probe_vel_kms,ads_accuracy_valid,precision_confirmed,'// &
            'first_failure_hour,first_unconfirmed_hour'
        write(output%stats_unit,'(A)') header
    end subroutine

    subroutine probability_write_hour(output,population,hour,epoch,states,reference,owners,valid,confirmed, &
        split_event,depth_event,failure_event,patches,depth_limited,original_pos,original_vel,probe_pos,probe_vel,is_final,status)
        type(probability_output_type),intent(inout) :: output
        type(probability_population_type),intent(in) :: population
        real(DP),intent(in) :: hour,epoch,states(:,:),reference(6),original_pos,original_vel,probe_pos,probe_vel
        integer,intent(in) :: owners(:),patches,depth_limited
        logical,intent(in) :: valid,confirmed,split_event,depth_event,failure_event,is_final
        integer,intent(out) :: status
        type(srp_mc_group_stats_type) :: stats
        integer :: n,p,i,j,unit
        logical :: skew_valid(6),kurt_valid(6)
        real(DP) :: moment_scale(6),population_sigma(6)
        character(len=256) :: kind
        status=-1; n=size(population%q,2)
        if(any(shape(states)/=[6,n]).or.size(owners)/=n) return
        if(.not.all(ieee_is_finite(states)).or..not.all(ieee_is_finite(reference)).or.any(owners<1)) return
        if(.not.all(ieee_is_finite([hour,epoch,original_pos,original_vel,probe_pos,probe_vel]))) return
        if(hour<=output%previous_hour) return
        if(.not.valid.and.output%first_failure_hour<0) output%first_failure_hour=hour
        if(.not.confirmed.and.output%first_unconfirmed_hour<0) output%first_unconfirmed_hour=hour
        unit=output%stats_unit
        do p=1,size(output%prefix_counts)
            call compute_group_stats(states(:,1:output%prefix_counts(p)),reference,stats,status)
            if(status/=0) return
            moment_scale=max(1.0_DP,abs(stats%component_min),abs(stats%component_max),abs(stats%mean_error))
            population_sigma=stats%standard_deviation*sqrt(real(stats%sample_count-1,DP)/real(stats%sample_count,DP))
            skew_valid=population_sigma>100*epsilon(1.0_DP)*moment_scale.and.stats%sample_count>2
            kurt_valid=population_sigma>100*epsilon(1.0_DP)*moment_scale.and.stats%sample_count>3
            write(unit,'(ES25.16E3)',advance='no') hour
            call values(unit,[epoch]); call integer_value(unit,stats%sample_count)
            call integer_value(unit,size(output%visual_ids))
            call values(unit,reference); call values(unit,reference+stats%mean_error); call values(unit,stats%mean_error)
            do i=1,6
                call values(unit,stats%covariance(i,:))
            end do
            call values(unit,stats%standard_deviation); call values(unit,stats%component_min)
            call values(unit,stats%component_max)
            do i=1,6
                call optional_value(unit,stats%skewness(i),skew_valid(i))
            end do
            do i=1,6
                call optional_value(unit,stats%excess_kurtosis(i),kurt_valid(i))
            end do
            do i=1,6
                call integer_value(unit,merge(1,0,skew_valid(i)))
            end do
            do i=1,6
                call integer_value(unit,merge(1,0,kurt_valid(i)))
            end do
            call values(unit,stats%position_axes_rms); call values(unit,stats%velocity_axes_rms)
            call values(unit,[stats%position_norm_rms,stats%velocity_norm_rms])
            call values(unit,stats%position_norm_quantiles); call values(unit,[stats%position_norm_max])
            call integer_value(unit,stats%position_worst_id)
            call values(unit,stats%velocity_norm_quantiles); call values(unit,[stats%velocity_norm_max])
            call integer_value(unit,stats%velocity_worst_id)
            call integer_value(unit,patches); call integer_value(unit,depth_limited)
            call values(unit,[original_pos,original_vel,probe_pos,probe_vel])
            call integer_value(unit,merge(1,0,valid)); call integer_value(unit,merge(1,0,confirmed))
            call optional_value(unit,output%first_failure_hour,output%first_failure_hour>=0)
            call optional_value(unit,output%first_unconfirmed_hour,output%first_unconfirmed_hour>=0)
            write(unit,*)
        end do
        call write_points(output%visual_unit,hour,epoch,'hourly_visualization',states,reference,owners,output%visual_ids)
        kind=''
        if(output%previous_hour<0) call label(kind,'start')
        if(is_final) call label(kind,'end')
        if(any(abs(hour-[0.0_DP,24.0_DP,72.0_DP,168.0_DP,240.0_DP,360.0_DP])<1e-10_DP)) call label(kind,'fixed')
        if(split_event.and..not.output%seen_split) then
            call label(kind,'first_split'); call label(output%previous_kind,'before_first_split')
            output%seen_split=.true.
        end if
        if(depth_event.and..not.output%seen_depth) then
            call label(kind,'first_depth_limit'); call label(output%previous_kind,'before_first_depth_limit')
            output%seen_depth=.true.
        end if
        if(failure_event.and..not.output%seen_failure) then
            call label(kind,'first_validation_failure'); call label(output%previous_kind,'before_first_validation_failure')
            output%seen_failure=.true.
        end if
        call flush_previous(output)
        output%previous_hour=hour; output%previous_epoch=epoch; output%previous_reference=reference
        output%previous_states=states; output%previous_owners=owners; output%previous_kind=kind
        status=0
    end subroutine

    subroutine probability_close(output,status)
        type(probability_output_type),intent(inout) :: output
        integer,intent(out) :: status
        integer :: io
        status=0
        if(output%key_unit/=0) call flush_previous(output)
        if(output%stats_unit/=0) then
            close(output%stats_unit,iostat=io); if(io/=0) status=io
        end if
        if(output%visual_unit/=0) then
            close(output%visual_unit,iostat=io); if(io/=0) status=io
        end if
        if(output%key_unit/=0) then
            close(output%key_unit,iostat=io); if(io/=0) status=io
        end if
        output%stats_unit=0; output%visual_unit=0; output%key_unit=0
    end subroutine

    subroutine flush_previous(output)
        type(probability_output_type),intent(inout) :: output
        if(output%previous_hour<0.or.len_trim(output%previous_kind)==0) return
        call write_points(output%key_unit,output%previous_hour,output%previous_epoch,trim(output%previous_kind), &
            output%previous_states,output%previous_reference,output%previous_owners,output%all_ids)
        output%previous_kind=''
    end subroutine

    subroutine write_points(unit,hour,epoch,kind,states,reference,owners,ids)
        integer,intent(in) :: unit,owners(:),ids(:)
        real(DP),intent(in) :: hour,epoch,states(:,:),reference(6)
        character(len=*),intent(in) :: kind
        integer :: i,j
        do j=1,size(ids)
            i=ids(j)
            write(unit,'(2(ES25.16E3,","),A,2(",",I0),12(",",ES25.16E3))') &
                hour,epoch,kind,i,owners(i),states(:,i),states(:,i)-reference
        end do
    end subroutine

    subroutine values(unit,array)
        integer,intent(in) :: unit
        real(DP),intent(in) :: array(:)
        integer :: i
        do i=1,size(array)
            write(unit,'(",",ES25.16E3)',advance='no') array(i)
        end do
    end subroutine
    subroutine integer_value(unit,value)
        integer,intent(in) :: unit,value
        write(unit,'(",",I0)',advance='no') value
    end subroutine
    subroutine optional_value(unit,value,valid)
        integer,intent(in) :: unit
        real(DP),intent(in) :: value
        logical,intent(in) :: valid
        if(valid) then
            call values(unit,[value])
        else
            write(unit,'(A)',advance='no') ','
        end if
    end subroutine
    subroutine label(text,item)
        character(len=*),intent(inout) :: text
        character(len=*),intent(in) :: item
        if(len_trim(text)==0) then
            text=item
        else
            text=trim(text)//'|'//item
        end if
    end subroutine
    function digit(i) result(text)
        integer,intent(in) :: i
        character(len=1) :: text
        write(text,'(I1)') i
    end function
    subroutine header_vector(header,name,n)
        character(len=:),allocatable,intent(inout) :: header
        character(len=*),intent(in) :: name
        integer,intent(in) :: n
        integer :: i
        do i=1,n
            header=header//','//name//digit(i)
        end do
    end subroutine
end module
