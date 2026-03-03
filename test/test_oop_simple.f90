program test_oop_simple
    use pod_global, only: DP
    use pod_object_base, only: pod_object
    use pod_system_oop_simple, only: pod_system_simple, pod_logger_simple
    implicit none
    
    ! Object instances
    type(pod_system_simple) :: system
    type(pod_logger_simple) :: logger
    
    write(*, *) '=========================================='
    write(*, *) '     CAT OOP Simple Test Program'
    write(*, *) '=========================================='
    write(*, *)
    
    ! Test 1: Logger functionality
    write(*, *) 'Test 1: Logger Functionality'
    call logger%set_name('Test_Logger_Simple')
    call logger%initialize()
    
    if (logger%is_initialized()) then
        write(*, *) 'Logger initialized successfully'
        write(*, *) 'Logger name: ', trim(logger%get_name())
        write(*, *) 'Log level: ', logger%get_log_level()
        write(*, *) 'Log file: ', trim(logger%get_log_file())
        write(*, *) 'Console output: ', logger%get_console_output()
        write(*, *) 'File output: ', logger%get_file_output()
        
        ! Test logging
        call logger%log_message(2, 'This is an info message', 'test_program')
        call logger%log_message(3, 'This is a warning message', 'test_program')
        call logger%log_message(4, 'This is an error message', 'test_program')
        
        ! Test log level change
        call logger%set_log_level(3)
        write(*, *) 'Log level changed to: ', logger%get_log_level()
        
        call logger%cleanup()
        write(*, *) 'Logger cleaned up'
    else
        write(*, *) 'Failed to initialize logger'
    end if
    write(*, *)
    
    ! Test 2: Main System
    write(*, *) 'Test 2: Main System'
    call system%set_name('Test_POD_System_Simple')
    call system%initialize()
    
    if (system%is_initialized()) then
        write(*, *) 'System initialized successfully'
        write(*, *) 'System name: ', trim(system%get_name())
        write(*, *) 'System state: ', trim(system%get_system_state_string())
        write(*, *) 'Memory usage: ', system%get_memory_usage(), ' MB'
        write(*, *) 'Error count: ', system%get_error_count()
        write(*, *) 'SPICE loaded: ', system%is_spice_loaded()
        
        ! Test system status
        call system%print_system_status()
        
        ! Test SPICE operations
        call system%load_spice_kernels()
        write(*, *) 'SPICE loaded: ', system%is_spice_loaded()
        
        call system%unload_spice_kernels()
        write(*, *) 'SPICE unloaded: ', system%is_spice_loaded()
        
        ! Test error handling
        if (system%has_error()) then
            write(*, *) 'System has error: ', trim(system%get_error_message())
        else
            write(*, *) 'System has no errors'
        end if
        
        call system%cleanup()
        write(*, *) 'System cleaned up'
    else
        write(*, *) 'Failed to initialize system'
    end if
    write(*, *)
    
    ! Test 3: Polymorphism demonstration
    write(*, *) 'Test 3: Polymorphism Demonstration'
    call test_polymorphism()
    write(*, *)
    
    write(*, *) '=========================================='
    write(*, *) '     All OOP Simple Tests Completed'
    write(*, *) '=========================================='
    
contains
    
    subroutine test_polymorphism()
        implicit none
        
        class(pod_object), allocatable :: obj1, obj2
        
        ! Allocate different types
        allocate(pod_system_simple::obj1)
        allocate(pod_logger_simple::obj2)
        
        ! Set names
        call obj1%set_name('Polymorphic_System_Simple')
        call obj2%set_name('Polymorphic_Logger_Simple')
        
        ! Test polymorphic behavior
        write(*, *) 'Polymorphic object 1 name: ', trim(obj1%get_name())
        write(*, *) 'Polymorphic object 2 name: ', trim(obj2%get_name())
        
        ! Test initialization
        call obj1%initialize()
        call obj2%initialize()
        
        write(*, *) 'Object 1 initialized: ', obj1%is_initialized()
        write(*, *) 'Object 2 initialized: ', obj2%is_initialized()
        
        ! Cleanup
        call obj1%cleanup()
        call obj2%cleanup()
        
        deallocate(obj1)
        deallocate(obj2)
        
        write(*, *) 'Polymorphism test completed'
    end subroutine test_polymorphism

end program test_oop_simple
