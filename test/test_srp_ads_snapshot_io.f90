program test_srp_ads_snapshot_io
    use pod_global, only: DP
    use pod_uq_ads_coordinates_module, only: ADS_COORD_WHITENED
    use pod_uq_hfem_ads_module, only: hfem_ads_options_type, &
        hfem_ads_history_type, hfem_ads_history_init, &
        hfem_ads_history_evaluate, hfem_ads_history_destroy
    use pod_srp_ads_snapshot_module, only: write_ads_patch_snapshot
    implicit none
    character(len=*), parameter :: prefix='/tmp/pod_srp_ads_snapshot_test'
    type(hfem_ads_options_type) :: options
    type(hfem_ads_history_type) :: history
    real(DP) :: nominal(6), covariance(6,6), point(7,1), expected(6,1), rebuilt(6)
    real(DP) :: coefficient, monomial
    integer :: status, unit, io, patch_id, component, powers(7), i
    logical :: found(1), saw_constant, saw_linear
    character(len=256) :: message, line

    nominal = [10.0_DP,20.0_DP,30.0_DP,1.0_DP,2.0_DP,3.0_DP]
    covariance = 0.0_DP
    do i=1,6
        covariance(i,i)=0.01_DP
    end do
    options%coordinate_mode=ADS_COORD_WHITENED
    options%srp_sigma=0.02_DP
    call hfem_ads_history_init(history,nominal,covariance,0.0_DP,options,status,message)
    call require(status==0,'initialize snapshot domain')
    point(:,1)=[0.2_DP,-0.3_DP,0.1_DP,0.0_DP,0.4_DP,-0.2_DP,0.5_DP]
    call hfem_ads_history_evaluate(history,point,expected,found,status)
    call require(status==0.and.all(found),'evaluate source patch')

    call write_ads_patch_snapshot(history,prefix,.false.,status,message)
    call require(status==0,'export snapshot: '//trim(message))
    rebuilt=0.0_DP
    saw_constant=.false.
    saw_linear=.false.
    open(newunit=unit,file=prefix//'_coeff.csv',status='old',action='read',iostat=io)
    call require(io==0,'open exported coefficients')
    read(unit,'(A)',iostat=io) line
    do
        read(unit,'(A)',iostat=io) line
        if(io/=0) exit
        read(line,*,iostat=io) patch_id,component,powers,coefficient
        call require(io==0,'parse sparse coefficient row')
        monomial=coefficient
        do i=1,7
            monomial=monomial*point(i,1)**powers(i)
        end do
        rebuilt(component)=rebuilt(component)+monomial
        if(component==1.and.all(powers==0)) saw_constant=.true.
        if(component==1.and.powers(1)==1.and.sum(powers)==1) saw_linear=.true.
    end do
    close(unit)
    call require(saw_constant.and.saw_linear,'constant and first-order terms exported')
    call require(maxval(abs(rebuilt-expected(:,1)))<1.0e-10_DP, &
        'exported polynomial reconstructs patch values')
    call write_ads_patch_snapshot(history,prefix,.true.,status,message)
    call require(status/=0,'previous snapshot unavailable at initial time')
    call hfem_ads_history_destroy(history)
    write(*,'(A)') 'PASS: sparse ADS snapshot reconstructs the initial patch'
contains
    subroutine require(ok,label)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: label
        if(ok) return
        write(*,'(A)') 'FAIL: '//trim(label)
        stop 1
    end subroutine require
end program test_srp_ads_snapshot_io
