!> @file test_compiled_da_batch_eval.f90
!! @brief Specify one-call, allocation-reusing evaluation of many DA points.
program test_compiled_da_batch_eval
    use pod_global,only: DP
    use pod_dace_classes,only: DA,AlgebraicVector,CompiledDA,dace_initialize, &
        active_da_count,da_add,da_sub,da_mul
    implicit none

    integer,parameter :: N_POINTS=4096
    type(DA) :: x1,x2,x3,t1,t2,t3
    type(AlgebraicVector) :: map
    type(CompiledDA) :: compiled,invalid
    real(DP) :: inputs(3,N_POINTS),batch_results(2,N_POINTS)
    real(DP) :: scalar_result(2),wrong_rows(1,N_POINTS)
    real(DP) :: wrong_points(2,N_POINTS-1),empty_inputs(3,0),empty_results(2,0)
    real(DP) :: max_difference
    integer :: point,status,handles_before

    call dace_initialize(2,3)
    call x1%init_var(1)
    call x2%init_var(2)
    call x3%init_var(3)
    call t1%init(); call t2%init(); call t3%init()
    call map%init(2)

    ! map(1)=1+x1+2*x2+x3^2
    call da_mul(x3,x3,t1)
    call da_mul(x2,2.0_DP,t2)
    call da_add(x1,t2,t3)
    call da_add(t3,t1,t2)
    call da_add(t2,1.0_DP,map%elements(1))
    ! map(2)=x1*x2-x3
    call da_mul(x1,x2,t1)
    call da_sub(t1,x3,map%elements(2))

    call x1%destroy(); call x2%destroy(); call x3%destroy()
    call t1%destroy(); call t2%destroy(); call t3%destroy()
    compiled=map%compile()
    do point=1,N_POINTS
        inputs(1,point)=-1.0_DP+2.0_DP*real(point-1,DP)/real(N_POINTS-1,DP)
        inputs(2,point)=sin(0.01_DP*real(point,DP))
        inputs(3,point)=cos(0.02_DP*real(point,DP))
    end do

    handles_before=active_da_count()
    call compiled%eval_batch_into(inputs,batch_results,status)
    call assert_true(status==0,'valid batch status')
    call assert_true(active_da_count()==handles_before, &
                     'batch evaluation must not create DA handles')
    max_difference=0.0_DP
    do point=1,N_POINTS
        call compiled%eval_into(inputs(:,point),scalar_result,status)
        call assert_true(status==0,'scalar reference status')
        max_difference=max(max_difference, &
                           maxval(abs(batch_results(:,point)-scalar_result)))
    end do
    call assert_close(max_difference,0.0_DP,1.0e-14_DP, &
                      'batch results equal scalar results')

    call compiled%eval_batch_into(inputs,wrong_rows,status)
    call assert_true(status==-1,'wrong result row count is rejected')
    call compiled%eval_batch_into(inputs,wrong_points,status)
    call assert_true(status==-1,'mismatched point count is rejected')
    call compiled%eval_batch_into(empty_inputs,empty_results,status)
    call assert_true(status==-1,'empty batch is rejected')
    call invalid%eval_batch_into(inputs,batch_results,status)
    call assert_true(status==-2,'invalid compiled handle is rejected')

    call compiled%destroy()
    call map%destroy()
    call assert_true(active_da_count()==0,'all DA handles are released')
    write(*,'(a)') 'PASS: CompiledDA batch evaluation matches scalar evaluation.'

contains

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
            write(*,'(a,3(1x,es24.16))') 'FAIL: '//trim(label), &
                actual,expected,tolerance
            error stop 1
        end if
    end subroutine assert_close

end program test_compiled_da_batch_eval
