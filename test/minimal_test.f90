program minimal_test
    use pod_global
    implicit none
    
    write(*, *) '最小测试 CAT Global 模块'
    
    ! 只测试初始化
    call pod_init()
    call pod_cleanup()
    
    write(*, *) '测试完成'
    
end program minimal_test
