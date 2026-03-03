!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!> # CAT SPICE Interface Module
!> 
!> This module provides a minimal Fortran interface to the NASA SPICE toolkit
!> for the CAT Fortran space object monitoring system. It focuses on essential
!> SPICE functionality for kernel management and basic interface declarations.
!> 
!> ## Features
!> 
!> - **SPICE Initialization**: Automatic loading of required SPICE kernels
!> - **Kernel Management**: Loading and cleanup of SPICE kernels
!> - **Interface Declarations**: Fortran interfaces to essential SPICE functions
!> 
!> ## SPICE Kernels Used
!> 
!> - **Leap Second Kernel**: `naif0012.tls` - For time conversions
!> - **Planetary Constants**: `pck00010.tpc` - For planetary parameters
!> - **Planetary Ephemeris**: `de421.bsp` - For planetary positions
!> - **Earth Orientation**: `earth_000101_230801_230509.bpc` - For Earth orientation
!> 
!> ## Dependencies
!> 
!> - **External**: NASA SPICE Toolkit (`spicelib.a`)
!> - **Internal**: `pod_global` for file operations and constants
!> 
!> ## Input/Output Files
!> 
!> - **Input**: SPICE kernel files in `kernels/` directory
!> - **Output**: Console messages for initialization and warnings
!> 
!> ## Author
!> 
!> Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> 
!> ## Version
!> 
!> - **Created**: 2025-09-08
!> - **Updated**: 2025-09-10
!> - **Ref**: 1. NASA SPICE Toolkit Documentation: https://naif.jpl.nasa.gov/naif/toolkit.html
!>            2. SPICE User's Guide: https://naif.jpl.nasa.gov/naif/tutorials.html
!>            3. SPICE Reference Manual: https://naif.jpl.nasa.gov/naif/reference_manual.html
!> 
!> @note This module serves as the primary interface between CAT Fortran and the SPICE toolkit.
!>       It provides only essential kernel management and interface declarations.
!>       Other modules can directly call SPICE functions using the declared interfaces.
!> 
!> @warning SPICE kernels must be available in the `kernels/` directory for proper operation.
!>          Missing kernels will result in warnings but will not prevent module initialization.
!> 
!> @todo Add support for custom kernel loading
!> @todo Add support for kernel validation
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------

module pod_spice
    use pod_global, only: DP, MAX_STRING_LEN, file_exists
    implicit none
    
    ! SPICE状态
    logical :: spice_initialized = .false.
    character(len=MAX_STRING_LEN) :: spice_kernel_path = './kernels/'
    
    ! SPICE常量
    integer, parameter :: SPICE_MAX_STRING_LEN = 256
    integer, parameter :: SPICE_MAX_KERNELS = 100
    
    ! SPICE函数接口声明
    interface
        subroutine furnsh(kernel)
            character(len=*), intent(in) :: kernel
        end subroutine furnsh
        
        subroutine unload(kernel)
            character(len=*), intent(in) :: kernel
        end subroutine unload
        
        subroutine spkezr(target, et, ref, abcorr, obs, state, found)
            character(len=*), intent(in) :: target
            real(8), intent(in) :: et
            character(len=*), intent(in) :: ref
            character(len=*), intent(in) :: abcorr
            character(len=*), intent(in) :: obs
            real(8), dimension(6), intent(out) :: state
            logical, intent(out) :: found
        end subroutine spkezr
        
        subroutine str2et(str, et)
            character(len=*), intent(in) :: str
            real(8), intent(out) :: et
        end subroutine str2et
        
        subroutine et2utc(et, format, prec, utcstr)
            real(8), intent(in) :: et
            character(len=*), intent(in) :: format
            integer, intent(in) :: prec
            character(len=*), intent(out) :: utcstr
        end subroutine et2utc
        
        subroutine pxform(from, to, et, rot)
            character(len=*), intent(in) :: from
            character(len=*), intent(in) :: to
            real(8), intent(in) :: et
            real(8), dimension(3,3), intent(out) :: rot
        end subroutine pxform
        
        subroutine bodvrd(body, item, maxn, dim, values)
            character(len=*), intent(in) :: body
            character(len=*), intent(in) :: item
            integer, intent(in) :: maxn
            integer, intent(out) :: dim
            double precision, intent(out) :: values(*)
        end subroutine bodvrd
        
        subroutine bodn2c(name, code, found)
            character(len=*), intent(in) :: name
            integer, intent(out) :: code
            logical, intent(out) :: found
        end subroutine bodn2c
        
        subroutine bodc2n(code, name, found)
            integer, intent(in) :: code
            character(len=*), intent(out) :: name
            logical, intent(out) :: found
        end subroutine bodc2n
        
        subroutine ttrans(input_type, output_type, time)
            character(len=*), intent(in) :: input_type
            character(len=*), intent(in) :: output_type
            real(8), intent(inout) :: time
        end subroutine ttrans
    end interface
    
contains

    !> Initialize the SPICE system by loading required kernels
    !> 
    !> This subroutine automatically loads the essential SPICE kernels needed for
    !> astronomical computations. It loads leap second, planetary constants,
    !> planetary ephemeris, and Earth orientation kernels.
    !> 
    !> @note This subroutine is called automatically by other SPICE functions
    !>       if the system is not already initialized.
    !> 
    !> @warning Missing kernel files will result in warnings but will not prevent
    !>          initialization from completing.
    subroutine spice_init()
        implicit none
        character(len=MAX_STRING_LEN) :: kernel_file
        
        if (spice_initialized) return
        
        write(*, *) '初始化SPICE系统...'
        
        ! 加载闰秒内核
        kernel_file = trim(spice_kernel_path) // 'lsk/naif0012.tls'
        if (file_exists(kernel_file)) then
            call furnsh(kernel_file)
            write(*, *) '已加载闰秒内核: ', trim(kernel_file)
        else
            write(*, *) '警告: 闰秒内核文件未找到: ', trim(kernel_file)
        end if
        
        ! 加载行星内核
        kernel_file = trim(spice_kernel_path) // 'pck/pck00010.tpc'
        if (file_exists(kernel_file)) then
            call furnsh(kernel_file)
            write(*, *) '已加载行星内核: ', trim(kernel_file)
        else
            write(*, *) '警告: 行星内核文件未找到: ', trim(kernel_file)
        end if
        
        ! 加载星历内核
        kernel_file = trim(spice_kernel_path) // 'spk/de421.bsp'
        if (file_exists(kernel_file)) then
            call furnsh(kernel_file)
            write(*, *) '已加载星历内核: ', trim(kernel_file)
        else
            write(*, *) '警告: 星历内核文件未找到: ', trim(kernel_file)
        end if
        
        ! 加载地球内核
        kernel_file = trim(spice_kernel_path) // 'pck/earth_000101_230801_230509.bpc'
        if (file_exists(kernel_file)) then
            call furnsh(kernel_file)
            write(*, *) '已加载地球内核: ', trim(kernel_file)
        else
            write(*, *) '警告: 地球内核文件未找到: ', trim(kernel_file)
        end if
        
        spice_initialized = .true.
        write(*, *) 'SPICE系统初始化完成'
    end subroutine spice_init
    
    !> Clean up SPICE resources and unload all kernels
    !> 
    !> This subroutine unloads all SPICE kernels and resets the initialization
    !> status. It should be called when the SPICE system is no longer needed.
    !> 
    !> @note This subroutine is safe to call even if SPICE is not initialized.
    subroutine spice_cleanup()
        implicit none
        
        if (.not. spice_initialized) return
        
        write(*, *) '清理SPICE资源...'
        call unload('*')
        spice_initialized = .false.
        write(*, *) 'SPICE资源清理完成'
    end subroutine spice_cleanup
    
    !> Check if SPICE system is initialized
    !> 
    !> @return .true. if SPICE is initialized, .false. otherwise
    logical function is_spice_initialized()
        is_spice_initialized = spice_initialized
    end function is_spice_initialized
    
    !> Set SPICE kernel path
    !> 
    !> @param[in] path New kernel path
    subroutine set_spice_kernel_path(path)
        character(len=*), intent(in) :: path
        spice_kernel_path = trim(path)
    end subroutine set_spice_kernel_path
    
    !> Get SPICE kernel path
    !> 
    !> @return Current kernel path
    function get_spice_kernel_path() result(path)
        character(len=MAX_STRING_LEN) :: path
        path = spice_kernel_path
    end function get_spice_kernel_path

end module pod_spice