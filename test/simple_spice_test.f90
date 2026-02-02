program simple_spice_test
    implicit none
    
    ! SPICE函数接口声明
    interface
        subroutine bodvrd(body, item, maxn, dim, values)
            character(len=*), intent(in) :: body
            character(len=*), intent(in) :: item
            integer, intent(in) :: maxn
            integer, intent(out) :: dim
            double precision, intent(out) :: values(*)
        end subroutine bodvrd
        
        subroutine furnsh(file)
            character(len=*), intent(in) :: file
        end subroutine furnsh
        
        subroutine unload(file)
            character(len=*), intent(in) :: file
        end subroutine unload
    end interface
    
    double precision :: gm_values(1)
    integer :: dim
    logical :: found
    
    write(*, *) '=========================================='
    write(*, *) '     简单SPICE测试程序'
    write(*, *) '=========================================='
    write(*, *)
    
    ! 尝试加载SPICE内核
    write(*, *) '1. 尝试加载SPICE内核...'
    call furnsh('kernels/lsk/naif0012.tls')
    call furnsh('kernels/pck/gm_de431.tpc')
    write(*, *) '   SPICE内核加载成功'
    
    ! 尝试读取地球的GM值
    write(*, *) '2. 尝试读取地球的GM值...'
    call bodvrd('EARTH', 'GM', 1, dim, gm_values)
    write(*, *) '   地球GM值: ', gm_values(1), ' km³/s²'
    
    ! 清理SPICE内核
    write(*, *) '3. 清理SPICE内核...'
    call unload('kernels/lsk/naif0012.tls')
    call unload('kernels/pck/gm_de431.tpc')
    write(*, *) '   SPICE内核清理完成'
    
    write(*, *)
    write(*, *) '=========================================='
    write(*, *) '     测试完成'
    write(*, *) '=========================================='

end program simple_spice_test
