!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!> # CAT Global Module
!> 
!> This module provides global definitions, constants, and system management
!> functions for the CAT Fortran space object monitoring system.
!> 
!> ## Features
!> 
!> - **Physical Constants**: Gravitational parameters, astronomical units, etc.
!> - **Error Handling**: Comprehensive error handling and logging system
!> - **System State Management**: System initialization and monitoring
!> - **Logging System**: Multi-level logging with file and console output
!> 
!> ## Dependencies
!> 
!> - None (this is a base module)
!> 
!> ## Input/Output Files
!> 
!> - **Input**: `config/pod_config.txt` (optional)
!> - **Output**: `logs/cat.log`
!> 
!> ## Author
!> 
!> Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> 
!> ## Version
!> 
!> - **Created**: 2025-09-08
!> - **Updated**: 2025-09-10  update: get_file_size, directory_exists
!> - **Ref**: 1. JPL Astrodynamic Parameters Website: https://ssd.jpl.nasa.gov/astro_par.html
!>            2. SPICE Toolkit Documentation: https://naif.jpl.nasa.gov/naif/toolkit.html
!>            3. WGS84 Standard: World Geodetic System 1984
!>            4. CODATA 2018: NIST Special Publication 330
!>            5. IERS Conventions 2010: https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

!> 
!> @note This module serves as the foundation for the entire CAT Fortran system
!> and should be initialized before using any other modules.
!> 
!> @warning Do not modify global constants without understanding their impact
!> on the entire system.
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------


module pod_global
    implicit none
    
    ! Make essential constants and functions public
    public :: DP, SP, MAX_STRING_LEN, MAX_ARRAY_SIZE
    public :: SUCCESS, ERROR_INVALID_INPUT, ERROR_FILE_NOT_FOUND, ERROR_MEMORY_ALLOCATION
    public :: LOG_DEBUG, LOG_INFO, LOG_WARNING, LOG_ERROR, LOG_CRITICAL
    public :: SPEED_OF_LIGHT, AU, EARTH_RADIUS, EARTH_MU, MOON_MU, SUN_MU
    public :: GRAVITATIONAL_CONSTANT
    public :: J2000_EPOCH, SECONDS_PER_DAY, DAYS_PER_YEAR
    public :: PLANET_MU
    public :: pod_init, pod_cleanup
    
    !> Double precision kind parameter
    !! Provides 15 decimal digits of precision and exponent range up to 10^307
    integer, parameter :: DP = selected_real_kind(15, 307)
    
    !> Single precision kind parameter  
    !! Provides 6 decimal digits of precision and exponent range up to 10^37
    integer, parameter :: SP = selected_real_kind(6, 37)
    
    !> Maximum string length for character variables
    integer, parameter :: MAX_STRING_LEN = 256
    
    !> Maximum array size for dynamic arrays
    integer, parameter :: MAX_ARRAY_SIZE = 10000
    
    !> Speed of light in vacuum
    !! Units: km/s
    !! Source: CODATA 2018 (exact value by definition)
    !! Reference: NIST Special Publication 330
    real(DP), parameter :: SPEED_OF_LIGHT = 299792.458_DP ! JPL
    
    !> Astronomical Unit - mean distance from Earth to Sun
    !! Units: km
    real(DP), parameter :: AU = 149597870.7_DP ! JPL
    
    !> Earth's equatorial radius
    !! Units: km
    !! Source: SPICE Toolkit (NASA JPL) - WGS84
    !! Reference: SPICE geophysical.ker file, WGS84 standard
    real(DP), parameter :: EARTH_RADIUS = 6378.137_DP ! this parameter might be determined differently under different conditions by different standards
    real(DP), parameter :: EARTH_RADIUS_b = 6356.752314245_DP ! short radius from WGS84
    
    !! Newton's constant of gravitation
    !! Units: km³/kg/s²
    !! Source: CODATA 2018 (exact value by definition)
    !! Reference: NIST Special Publication 330
    real(DP), parameter :: GRAVITATIONAL_CONSTANT = 6.67430e-11_DP ! JPL

    !> Earth's gravitational parameter (GM)
    !! Units: km³/s²
    !! Source: SPICE Toolkit (NASA JPL) - WGS84
    !! Reference: JPL Website: https://ssd.jpl.nasa.gov/astro_par.html
    !! Note: WGS84 value = 398600.5 km³/s²
    real(DP), parameter :: EARTH_MU = 398600.435507_DP ! JPL  
    !real(DP), parameter :: EARTH_MU = 398600.4418_DP ! WGS84
    
    !> Moon's gravitational parameter (GM)
    !! Units: km³/s²
    !! Source: SPICE Toolkit (NASA JPL)
    !! Reference: JPL Website: https://ssd.jpl.nasa.gov/astro_par.html
    !! Note: Need to confirm accurate values in SPICE
    real(DP), parameter :: MOON_MU = 4902.800118_DP ! JPL
    
    !> Sun's gravitational parameter (GM)
    !! Units: km³/s²
    !! Source: SPICE Toolkit (NASA JPL)
    !! Reference: JPL Website: https://ssd.jpl.nasa.gov/astro_par.html
    !! Note: Need to confirm accurate values in SPICE
    real(DP), parameter :: SUN_MU = 132712440041.279419_DP ! JPL

    !> Planetary gravitational parameters (GM)
    !! Units: km³/s²
    !! Source: SPICE Toolkit (NASA JPL)
    !! Reference: JPL Website: https://ssd.jpl.nasa.gov/astro_par.html
    !! Note: Need to confirm accurate values in SPICE
    !! Array: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Moon, Sun
    real(DP), parameter :: PLANET_MU(11) = (/ &
        22031.868551_DP, 324858.592000_DP, 398600.435507_DP, 42828.375816_DP, &
        126712764.1_DP, 37940584.8418_DP, 5794556.4_DP, 6836527.10058_DP, &
        975.5_DP, 4902.800118_DP, 132712440041.279419_DP /) ! JPL
    


    ! Time constants 
    real(DP), parameter :: SECONDS_PER_DAY = 86400.0_DP
    real(DP), parameter :: SECONDS_PER_HOUR = 3600.0_DP
    real(DP), parameter :: SECONDS_PER_MINUTE = 60.0_DP
    real(DP), parameter :: DAYS_PER_YEAR = 365.25_DP
    real(DP), parameter :: J2000_EPOCH = 2451545.0_DP  ! J2000.0 Julian Date


    ! PI and other mathematical constants
    real(DP), parameter :: PI = 3.1415926535897932_DP
    real(DP), parameter :: DEG_TO_RAD = PI / 180.0_DP
    real(DP), parameter :: RAD_TO_DEG = 180.0_DP / PI
    real(DP), parameter :: TWO_PI = 2.0_DP * PI
    real(DP), parameter :: HALF_PI = PI / 2.0_DP

    
    ! Error code definitions
    integer, parameter :: SUCCESS = 0
    integer, parameter :: ERROR_INVALID_INPUT = 1001
    integer, parameter :: ERROR_FILE_NOT_FOUND = 1002
    integer, parameter :: ERROR_FILE_IO = 1003
    integer, parameter :: ERROR_SPICE_ERROR = 1004
    integer, parameter :: ERROR_NUMERICAL_ERROR = 1005
    integer, parameter :: ERROR_MEMORY_ALLOCATION = 1006
    integer, parameter :: ERROR_CONVERGENCE = 1007
    integer, parameter :: ERROR_INVALID_CONFIG = 1008
    integer, parameter :: ERROR_SYSTEM_NOT_INITIALIZED = 1009
    
    ! Log levels
    integer, parameter :: LOG_DEBUG = 1
    integer, parameter :: LOG_INFO = 2
    integer, parameter :: LOG_WARNING = 3
    integer, parameter :: LOG_ERROR = 4
    integer, parameter :: LOG_CRITICAL = 5
    
    ! Global variables
    logical :: system_initialized = .false.
    character(len=MAX_STRING_LEN) :: spice_kernel_path
    character(len=MAX_STRING_LEN) :: output_directory
    character(len=MAX_STRING_LEN) :: log_file_path
    
    ! System state
    type system_state
        logical :: initialized = .false.
        logical :: spice_loaded = .false.
        logical :: config_loaded = .false.
        integer :: error_code = 0
        character(len=MAX_STRING_LEN) :: error_message = ''
        integer :: error_count = 0
        real(DP) :: memory_usage = 0.0_DP
        character(len=MAX_STRING_LEN) :: last_error = ''
        integer :: log_level = LOG_INFO
        logical :: verbose_output = .true.
    end type system_state
    
    type(system_state) :: global_state
    
    ! Log file unit number
    integer, parameter :: LOG_UNIT = 99
    
contains

    !> Initialize the CAT Fortran system
    !!
    !! This subroutine initializes the global system state, sets up logging,
    !! and prepares the system for use. It should be called before using
    !! any other CAT modules.
    !!
    !! @note This subroutine is idempotent - calling it multiple times
    !!       will not cause issues.
    !!
    !! @warning Must be called before using any other CAT modules
    subroutine pod_init()
        implicit none
        
        ! Check if already initialized
        if (system_initialized) then
            call log_message(LOG_WARNING, 'System already initialized', 'pod_init')
            return
        end if
        
        call log_message(LOG_INFO, 'Initializing CAT Fortran system...', 'pod_init')
        
        ! Initialize system state
        global_state%initialized = .false.
        global_state%spice_loaded = .false.
        global_state%config_loaded = .false.
        global_state%error_code = SUCCESS
        global_state%error_message = ''
        global_state%error_count = 0
        global_state%memory_usage = 0.0_DP
        global_state%last_error = ''
        global_state%log_level = LOG_INFO
        global_state%verbose_output = .true.
        
        ! Set default values
        spice_kernel_path = './kernels/'
        output_directory = './output/'
        log_file_path = './logs/cat.log'
        
        ! Create necessary directories
        call create_directories()
        
        ! Initialize logging system
        call init_logging()
        
        ! Note: SPICE initialization will be handled by other modules when needed
        
        system_initialized = .true.
        global_state%initialized = .true.
        
        call log_message(LOG_INFO, 'System initialization completed', 'pod_init')
    end subroutine pod_init
    
    !> Cleanup the CAT Fortran system
    !!
    !! This subroutine cleans up the global system state, closes log files,
    !! and resets all system variables. It should be called when the
    !! application is shutting down.
    !!
    !! @note This subroutine is idempotent - calling it multiple times
    !!       will not cause issues.
    !!
    !! @warning Should be called before program termination
    subroutine pod_cleanup()
        implicit none
        
        if (.not. system_initialized) then
            call log_message(LOG_WARNING, 'System not initialized, no cleanup needed', 'pod_cleanup')
            return
        end if
        
        call log_message(LOG_INFO, 'Cleaning up CAT Fortran system...', 'pod_cleanup')
        
        ! Note: SPICE cleanup will be handled by other modules
        
        ! Close log file
        call cleanup_logging()
        
        ! Reset system state
        system_initialized = .false.
        global_state%initialized = .false.
        global_state%spice_loaded = .false.
        global_state%config_loaded = .false.
        
        call log_message(LOG_INFO, 'System cleanup completed', 'pod_cleanup')
    end subroutine pod_cleanup
    
    subroutine create_directories()
        implicit none
        character(len=MAX_STRING_LEN) :: cmd
        
        ! Create output directory
        cmd = 'mkdir -p ' // trim(output_directory)
        call system(cmd)
        
        ! Create log directory
        cmd = 'mkdir -p ./logs'
        call system(cmd)
        
        ! Create config directory
        cmd = 'mkdir -p ./config'
        call system(cmd)
        
        ! Create kernel directory
        cmd = 'mkdir -p ' // trim(spice_kernel_path)
        call system(cmd)
    end subroutine create_directories
    
    ! SPICE-related functions will be handled by pod_spice module
    
    logical function file_exists(filename)
        character(len=*), intent(in) :: filename
        inquire(file=filename, exist=file_exists)
    end function file_exists
    
    subroutine set_error(error_code, error_message, module_name)
        integer, intent(in) :: error_code
        character(len=*), intent(in) :: error_message
        character(len=*), intent(in), optional :: module_name
        
        global_state%error_code = error_code
        global_state%error_message = error_message
        global_state%error_count = global_state%error_count + 1
        global_state%last_error = error_message
        
        ! Log error message
        if (present(module_name)) then
            call log_message(LOG_ERROR, 'Error [' // trim(module_name) // ']: ' // trim(error_message), 'set_error')
        else
            call log_message(LOG_ERROR, 'Error: ' // trim(error_message), 'set_error')
        end if
    end subroutine set_error
    
    subroutine clear_error()
        global_state%error_code = SUCCESS
        global_state%error_message = ''
        global_state%last_error = ''
    end subroutine clear_error
    
    logical function has_error()
        has_error = (global_state%error_code /= SUCCESS)
    end function has_error
    
    integer function get_error_code()
        get_error_code = global_state%error_code
    end function get_error_code
    
    function get_error_message() result(message)
        character(len=MAX_STRING_LEN) :: message
        message = global_state%error_message
    end function get_error_message
    
    ! Logging system
    subroutine init_logging()
        implicit none
        integer :: ios
        
        ! Open log file
        open(unit=LOG_UNIT, file=log_file_path, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*, *) 'Warning: Cannot create log file: ', trim(log_file_path)
            return
        end if
        
        ! Write log header
        write(LOG_UNIT, '(A)') '=========================================='
        write(LOG_UNIT, '(A)') '    CAT Fortran System Log'
        write(LOG_UNIT, '(A)') '=========================================='
        write(LOG_UNIT, *)
    end subroutine init_logging
    
    subroutine cleanup_logging()
        implicit none
        close(LOG_UNIT)
    end subroutine cleanup_logging
    
    subroutine log_message(level, message, module_name)
        integer, intent(in) :: level
        character(len=*), intent(in) :: message
        character(len=*), intent(in), optional :: module_name
        character(len=32) :: level_str, timestamp
        character(len=MAX_STRING_LEN) :: log_entry
        
        ! Check log level
        if (level < global_state%log_level) return
        
        ! Get timestamp
        call get_timestamp(timestamp)
        
        ! Set log level string
        select case (level)
            case (LOG_DEBUG)
                level_str = 'DEBUG'
            case (LOG_INFO)
                level_str = 'INFO '
            case (LOG_WARNING)
                level_str = 'WARN '
            case (LOG_ERROR)
                level_str = 'ERROR'
            case (LOG_CRITICAL)
                level_str = 'CRIT '
            case default
                level_str = 'UNKNW'
        end select
        
        ! Build log entry
        if (present(module_name)) then
            log_entry = '[' // trim(timestamp) // '] [' // trim(level_str) // '] [' // &
                       trim(module_name) // '] ' // trim(message)
        else
            log_entry = '[' // trim(timestamp) // '] [' // trim(level_str) // '] ' // trim(message)
        end if
        
        ! Write to log file
        write(LOG_UNIT, '(A)') trim(log_entry)
        flush(LOG_UNIT)
        
        ! If verbose output is enabled, also output to console
        if (global_state%verbose_output .and. level >= LOG_INFO) then
            write(*, '(A)') trim(log_entry)
        end if
    end subroutine log_message
    
    subroutine set_log_level(level)
        integer, intent(in) :: level
        global_state%log_level = level
        call log_message(LOG_INFO, 'Log level set to: ' // get_level_string(level), 'set_log_level')
    end subroutine set_log_level
    
    function get_level_string(level) result(level_str)
        integer, intent(in) :: level
        character(len=16) :: level_str
        
        select case (level)
            case (LOG_DEBUG)
                level_str = 'DEBUG'
            case (LOG_INFO)
                level_str = 'INFO'
            case (LOG_WARNING)
                level_str = 'WARNING'
            case (LOG_ERROR)
                level_str = 'ERROR'
            case (LOG_CRITICAL)
                level_str = 'CRITICAL'
            case default
                level_str = 'UNKNOWN'
        end select
    end function get_level_string
    
    subroutine get_timestamp(timestamp)
        character(len=32), intent(out) :: timestamp
        character(len=8) :: date
        character(len=10) :: time
        character(len=5) :: zone
        
        call date_and_time(date, time, zone)
        timestamp = date(1:4) // '-' // date(5:6) // '-' // date(7:8) // ' ' // &
                   time(1:2) // ':' // time(3:4) // ':' // time(5:6)
    end subroutine get_timestamp
    
    ! System status monitoring
    subroutine get_system_status(status)
        type(system_state), intent(out) :: status
        status = global_state
    end subroutine get_system_status
    
    subroutine print_system_status()
        implicit none
        
        write(*, *) '=== System Status ==='
        write(*, *) 'Initialized: ', global_state%initialized
        write(*, *) 'SPICE Loaded: ', global_state%spice_loaded
        write(*, *) 'Config Loaded: ', global_state%config_loaded
        write(*, *) 'Error Code: ', global_state%error_code
        write(*, *) 'Error Count: ', global_state%error_count
        write(*, *) 'Log Level: ', get_level_string(global_state%log_level)
        write(*, *) 'Verbose Output: ', global_state%verbose_output
        if (len_trim(global_state%last_error) > 0) then
            write(*, *) 'Last Error: ', trim(global_state%last_error)
        end if
        write(*, *)
    end subroutine print_system_status
    
    subroutine set_verbose_output(verbose)
        logical, intent(in) :: verbose
        global_state%verbose_output = verbose
        call log_message(LOG_INFO, 'Verbose output set to: ' // merge('ON ', 'OFF', verbose), 'set_verbose_output')
    end subroutine set_verbose_output
    
    logical function is_system_initialized()
        is_system_initialized = system_initialized
    end function is_system_initialized
    
    !> Get Earth's gravitational parameter
    !!
    !! @return Earth's gravitational parameter (GM) in km³/s²
    function get_earth_mu() result(mu)
        real(DP) :: mu
        mu = EARTH_MU
    end function get_earth_mu
    
    !> Get Moon's gravitational parameter
    !!
    !! @return Moon's gravitational parameter (GM) in km³/s²
    function get_moon_mu() result(mu)
        real(DP) :: mu
        mu = MOON_MU
    end function get_moon_mu
    
    !> Get Sun's gravitational parameter
    !!
    !! @return Sun's gravitational parameter (GM) in km³/s²
    function get_sun_mu() result(mu)
        real(DP) :: mu
        mu = SUN_MU
    end function get_sun_mu
    
    !> Get Earth's equatorial radius
    !!
    !! @return Earth's equatorial radius in km
    function get_earth_radius() result(radius)
        real(DP) :: radius
        radius = EARTH_RADIUS
    end function get_earth_radius
    
    !> Get Astronomical Unit
    !!
    !! @return Astronomical Unit (mean Earth-Sun distance) in km
    function get_au() result(au_value)
        real(DP) :: au_value
        au_value = AU
    end function get_au
    
    
    !> Get file size in bytes
    !> 
    !> This function returns the size of a file in bytes.
    !> 
    !> @param[in] filename Name of the file
    !> @return File size in bytes, or -1 if file doesn't exist
    integer function get_file_size(filename)
        character(len=*), intent(in) :: filename
        integer :: file_unit, ios
        
        if (.not. file_exists(filename)) then
            get_file_size = -1
            return
        end if
        
        open(newunit=file_unit, file=filename, status='old', iostat=ios)
        if (ios == 0) then
            inquire(unit=file_unit, size=get_file_size)
            close(file_unit)
        else
            get_file_size = -1
        end if
    end function get_file_size
    
    !> Check if a directory exists
    !> 
    !> This function checks whether a specified directory exists.
    !> 
    !> @param[in] dirname Name of the directory to check
    !> @return .true. if directory exists, .false. otherwise
    logical function directory_exists(dirname)
        character(len=*), intent(in) :: dirname
        logical :: exists
        inquire(file=trim(dirname)//'/.', exist=exists)
        directory_exists = exists
    end function directory_exists

end module pod_global
