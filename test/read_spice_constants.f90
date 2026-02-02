program read_spice_constants
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
    
    double precision :: values(10)
    integer :: dim
    
    write(*, *) '=========================================='
    write(*, *) '    SPICE 行星常数读取程序'
    write(*, *) '=========================================='
    write(*, *)
    
    ! 加载SPICE内核文件
    write(*, *) '加载SPICE内核文件...'
    call furnsh('kernels/lsk/naif0012.tls')
    call furnsh('kernels/pck/gm_de431.tpc')
    write(*, *) '内核文件加载完成'
    write(*, *)
    
    ! 读取地球引力参数
    write(*, *) '=== 地球参数 ==='
    call bodvrd('EARTH', 'GM', 1, dim, values)
    write(*, '(A, F15.6, A)') '地球引力参数 (GM): ', values(1), ' km³/s²'
    
    write(*, *)
    
    ! 读取月球参数
    write(*, *) '=== 月球参数 ==='
    call bodvrd('MOON', 'GM', 1, dim, values)
    write(*, '(A, F15.6, A)') '月球引力参数 (GM): ', values(1), ' km³/s²'
    
    write(*, *)
    
    ! 读取太阳参数
    write(*, *) '=== 太阳参数 ==='
    call bodvrd('SUN', 'GM', 1, dim, values)
    write(*, '(A, F15.6, A)') '太阳引力参数 (GM): ', values(1), ' km³/s²'
    
    write(*, *)
    
    ! 清理SPICE内核
    write(*, *) '清理SPICE内核...'
    call unload('kernels/lsk/naif0012.tls')
    call unload('kernels/pck/gm_de431.tpc')
    write(*, *) '清理完成'
    
    write(*, *)
    write(*, *) '=========================================='
    write(*, *) '    程序执行完成'
    write(*, *) '=========================================='
    
end program read_spice_constants