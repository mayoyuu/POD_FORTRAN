program test_srp_mc_io
    use pod_global, only: DP
    use pod_srp_mc_io_module, only: SRP_MC_POINTS_SCHEMA, read_srp_mc_samples, &
        write_mc_point_header, write_mc_point_rows
    implicit none

    real(DP),allocatable :: z(:,:),u(:,:)
    real(DP) :: treatment(6,2),control(6,2)
    integer :: ids(2),status,unit_value
    character(len=256) :: message,line
    character(len=*),parameter :: VALID='/tmp/pod_srp_mc_valid.csv'
    character(len=*),parameter :: BAD='/tmp/pod_srp_mc_bad.csv'
    character(len=*),parameter :: POINTS='/tmp/pod_srp_mc_points.csv'

    call write_samples(VALID,0)
    call read_srp_mc_samples(VALID,2,z,u,status,message)
    call assert_true(status==0,'valid read: '//trim(message))
    call assert_true(all(shape(z)==[7,2]),'z shape')
    call assert_true(all(shape(u)==[7,2]),'u shape')
    call assert_close(z(7,2),-2.1_DP,1.0e-14_DP,'z value')
    call assert_close(u(7,2),-0.7_DP,1.0e-14_DP,'u value')

    call write_samples(BAD,6)
    call read_srp_mc_samples(BAD,2,z,u,status,message,0.2_DP)
    call assert_true(status==0,'configurable eta sigma read: '//trim(message))
    call read_srp_mc_samples(BAD,2,z,u,status,message,0.3_DP)
    call assert_true(status/=0,'eta sigma mismatch rejected')

    call write_samples(BAD,1)
    call read_srp_mc_samples(BAD,2,z,u,status,message)
    call assert_true(status/=0,'duplicate id rejected')
    call write_samples(BAD,2)
    call read_srp_mc_samples(BAD,2,z,u,status,message)
    call assert_true(status/=0,'out of domain rejected')
    call write_samples(BAD,3)
    call read_srp_mc_samples(BAD,2,z,u,status,message)
    call assert_true(status/=0,'nonfinite rejected')
    call write_samples(BAD,4)
    call read_srp_mc_samples(BAD,2,z,u,status,message)
    call assert_true(status/=0,'malformed row rejected')
    call write_samples(BAD,5)
    call read_srp_mc_samples(BAD,2,z,u,status,message)
    call assert_true(status/=0,'wrong header rejected')
    call read_srp_mc_samples(VALID,3,z,u,status,message)
    call assert_true(status/=0,'wrong expected count rejected')

    treatment=reshape([(real(status,DP),status=1,12)],[6,2])
    control=0.5_DP*treatment
    ids=[11,12]
    open(newunit=unit_value,file=POINTS,status='replace',action='write')
    call write_mc_point_header(unit_value,status,message)
    call assert_true(status==0,'point header')
    call write_mc_point_rows(unit_value,2.0_DP,'fixed',ids,treatment,control,status,message)
    call assert_true(status==0,'point rows')
    close(unit_value)
    open(newunit=unit_value,file=POINTS,status='old',action='read')
    read(unit_value,'(A)') line
    call assert_true(index(trim(line),'schema='//SRP_MC_POINTS_SCHEMA)==1,'schema header')
    read(unit_value,'(A)') line
    call assert_true(index(line,'hour,keyframe_kind,sample_id')==1,'column header')
    read(unit_value,'(A)') line
    call assert_true(index(line,',fixed,11,')>0,'first point id and kind')
    close(unit_value)

    write(*,'(a)') 'PASS: SRP MC sample and point I/O.'

contains

    subroutine write_samples(path,mode)
        character(len=*),intent(in) :: path
        integer,intent(in) :: mode
        integer :: output_unit
        character(len=512) :: second_row

        open(newunit=output_unit,file=path,status='replace',action='write')
        if(mode==5) then
            write(output_unit,'(a)') 'sample_id,z1,z2'
        else
            write(output_unit,'(a)') 'sample_id,z1,z2,z3,z4,z5,z6,z7,'// &
                'u1,u2,u3,u4,u5,u6,u7,eta_treatment'
        end if
        write(output_unit,'(a)') '1,0.3,-0.6,0.9,-1.2,1.5,-1.8,2.1,'// &
            '0.1,-0.2,0.3,-0.4,0.5,-0.6,0.7,0.63'
        second_row='2,-0.3,0.6,-0.9,1.2,-1.5,1.8,-2.1,'// &
            '-0.1,0.2,-0.3,0.4,-0.5,0.6,-0.7,-0.63'
        select case(mode)
        case(1)
            second_row(1:1)='1'
        case(2)
            second_row='2,-3.1,0.6,-0.9,1.2,-1.5,1.8,-2.1,'// &
                '-1.0333333333333333,0.2,-0.3,0.4,-0.5,0.6,-0.7,-0.63'
        case(3)
            second_row='2,NaN,0.6,-0.9,1.2,-1.5,1.8,-2.1,'// &
                'NaN,0.2,-0.3,0.4,-0.5,0.6,-0.7,-0.63'
        case(4)
            second_row='2,-0.3,broken'
        case(6)
            second_row='2,-0.3,0.6,-0.9,1.2,-1.5,1.8,-2.1,'// &
                '-0.1,0.2,-0.3,0.4,-0.5,0.6,-0.7,-0.42'
        end select
        if(mode==6) then
            rewind(output_unit)
            write(output_unit,'(a)') 'sample_id,z1,z2,z3,z4,z5,z6,z7,'// &
                'u1,u2,u3,u4,u5,u6,u7,eta_treatment'
            write(output_unit,'(a)') '1,0.3,-0.6,0.9,-1.2,1.5,-1.8,2.1,'// &
                '0.1,-0.2,0.3,-0.4,0.5,-0.6,0.7,0.42'
        end if
        write(output_unit,'(a)') trim(second_row)
        close(output_unit)
    end subroutine write_samples

    subroutine assert_true(condition,label)
        logical,intent(in) :: condition
        character(len=*),intent(in) :: label
        if(.not.condition) then
            write(*,'(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual,expected,tolerance,label)
        real(DP),intent(in) :: actual,expected,tolerance
        character(len=*),intent(in) :: label
        if(abs(actual-expected)>tolerance) then
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label),actual,expected,tolerance
            error stop 1
        end if
    end subroutine assert_close

end program test_srp_mc_io
