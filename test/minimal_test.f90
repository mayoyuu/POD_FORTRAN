program minimal_test
    use cat_global
    implicit none
    
    write(*, *) '最小测试 CAT Global 模块'
    
    ! 只测试初始化
    call cat_init()
    call cat_cleanup()
    
    write(*, *) '测试完成'
    
end program minimal_test
