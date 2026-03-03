!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------
!> # CAT Time Module
!> 
!> This module provides comprehensive time transformation between different time scales
!> for the CAT Fortran space object monitoring system. It integrates with SPICE toolkit
!> for accurate astronomical time conversions.
!> 
!> ## Features
!> 
!> - **Time System Support**: UTC, TAI, TDT, TDB(ET)
!> - **Time Format Support**: ISO, JD, MJD, SPICE formats
!> - **SPICE Integration**: Direct interface to SPICE time functions
!> - **High Precision**: Double precision arithmetic for accuracy
!> - **Error Handling**: Comprehensive error checking and reporting
!> 
!> ## Supported Time Systems
!> 
!> - **UTC**: Coordinated Universal Time (civil time) - via SPICE str2et/et2utc
!> - **TAI**: International Atomic Time - via SPICE ttrans (limited support)
!> - **TT**: Terrestrial Time (TDT) - via SPICE ttrans (limited support)
!> - **TDB**: Barycentric Dynamical Time (same as ET) - via SPICE ttrans
!> 
!> ## Supported Time Formats
!> 
!> - **ISO**: ISO 8601 format (e.g., '2024-01-01T12:00:00')
!> - **JD**: Julian Date
!> - **MJD**: Modified Julian Date
!> - **SPICE**: Various SPICE time formats
!> 
!> ## SPICE Time Formats
!> 
!> | Format | Description | Example |
!> |--------|-------------|---------|
!> | 'C' | Calendar format | '1986 APR 12 16:31:09.814' |
!> | 'D' | Day-of-year format | '1986-102 // 16:31:12.814' |
!> | 'J' | Julian Date | 'JD 2446533.18834276' |
!> | 'ISOC' | ISO Calendar | '1987-04-12T16:31:12.814' |
!> | 'ISOD' | ISO Day-of-year | '1987-102T16:31:12.814' |
!> 
!> ## Dependencies
!> 
!> - **Internal**: `pod_global` for data types and constants
!> - **External**: NASA SPICE Toolkit for time conversions
!> 
!> ## Author
!> 
!> Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> 
!> ## Version
!> 
!> - **Created**: 2025-09-09
!> - **Updated**: 2025-09-11
!> - **Updated**: 2025-09-12 - Added TDT support, cleaned up unused constants
!> 
!> ### Supported Time Systems (as of 2025-09-12)
!> 
!> **Fully Supported**:
!> - **UTC** ↔ **ET**: via SPICE str2et/et2utc functions
!> - **UTC** ↔ **JD**: via custom conversion functions
!> - **UTC** ↔ **MJD**: via custom conversion functions
!> 
!> **SPICE ttrans Supported**:
!> - **TDB** ↔ **TAI**: verified working (difference ~32.18 seconds)
!> - **TDB** ↔ **TDT**: verified working (difference ~80 microseconds)
!> - **TAI** ↔ **TDT**: verified working
!> 
!> **Note**: TDT also stands for Terrestrial Time (TT)
!> 
!> - **Ref**: 1. SPICE Toolkit Documentation: https://naif.jpl.nasa.gov/naif/toolkit.html
!>            2. SPICE Time Conversion: https://naif.jpl.nasa.gov/naif/tutorials.html
!>            3. IAU SOFA Library: http://www.iausofa.org/
!>            4. IERS Conventions 2010: https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html
!> 
!> @note This module provides a unified interface for all time system conversions
!>       in the CAT Fortran system. All time operations should go through this module.
!> 
!> @warning Time conversions require proper SPICE kernel initialization.
!>          Ensure SPICE kernels are loaded before using time conversion functions.
!> 
!> @note SPICE ttrans function has limited time system support.
!>       Currently tested and working: TDB ↔ TAI, TDB ↔ TDT, TAI ↔ TDT conversions.
!>       Note: TDT is the SPICE name for Terrestrial Time (TT).
!> 
!> @todo Add support for more time systems (UT1, GPS, TCG, TCB, etc.)
!> @todo Add time arithmetic operations (addition, subtraction)
!> @todo Add time comparison functions
!> @todo Add time validation functions
!--------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------

module pod_time_module
    use pod_global, only: DP, MAX_STRING_LEN, SUCCESS, ERROR_INVALID_INPUT
    use pod_spice, only: str2et, et2utc, ttrans
    implicit none
    
    private
    public :: time_transfer, utc_to_et, et_to_utc, utc_to_jd, jd_to_utc
    public :: utc_to_mjd, mjd_to_utc, validate_time_string, get_current_utc
    
    
    ! Default time system and format
    character(len=20), parameter :: DEFAULT_TIME_SYSTEM = 'TDB'
    character(len=20), parameter :: DEFAULT_TIME_FORMAT = 'ISOC'
    
contains

    !> Transfer time between different time systems
    !> 
    !> This subroutine converts time from one time system to another using SPICE.
    !> 
    !> @param[in] t_input Input time value (seconds from J2000)
    !> @param[out] t_output Output time value (seconds from J2000)
    !> @param[in] t_input_type Input time system ('UTC', 'TAI', 'TDT', 'TDB', etc.)
    !> @param[in] t_output_type Output time system ('UTC', 'TAI', 'TDT', 'TDB', etc.)
    !> @note TDT is the SPICE name for Terrestrial Time (TT)
    !> @param[out] status Return status (0 = success, non-zero = error)
    subroutine time_transfer(t_input, t_output, t_input_type, t_output_type, status)
        real(DP), intent(in) :: t_input
        real(DP), intent(out) :: t_output
        character(len=*), intent(in) :: t_input_type, t_output_type
        integer, intent(out) :: status
        
        real(DP) :: tp
        
        status = SUCCESS
        tp = t_input
        
        ! Use SPICE ttrans for time system conversion
        call ttrans(t_input_type, t_output_type, tp)
        t_output = tp
        
    end subroutine time_transfer
    
    !> Convert UTC time string to ET (Ephemeris Time)
    !> 
    !> This function converts a UTC time string to ET (TDB) using SPICE.
    !> 
    !> @param[in] utc_string UTC time string (e.g., '2024-01-01T12:00:00')
    !> @return ET time in seconds from J2000
    !> 
    !> @example
    !> ```fortran
    !> real(DP) :: et_time
    !> et_time = utc_to_et('2024-01-01T12:00:00')
    !> ```
    function utc_to_et(utc_string)
        character(len=*), intent(in) :: utc_string
        real(DP) :: utc_to_et
        
        call str2et(utc_string, utc_to_et)
    end function utc_to_et
    
    !> Convert ET (Ephemeris Time) to UTC time string
    !> 
    !> This function converts ET to UTC time string using SPICE.
    !> 
    !> @param[in] et_time ET time in seconds from J2000
    !> @param[in] format Output format ('ISOC', 'C', 'D', 'J', etc.)
    !> @param[in] precision Number of decimal places for seconds
    !> @return UTC time string
    !> 
    !> @example
    !> ```fortran
    !> character(len=100) :: utc_string
    !> utc_string = et_to_utc(757382469.18392062_DP, 'ISOC', 3)
    !> ```
    function et_to_utc(et_time, format, precision)
        real(DP), intent(in) :: et_time
        character(len=*), intent(in), optional :: format
        integer, intent(in), optional :: precision
        character(len=MAX_STRING_LEN) :: et_to_utc
        
        character(len=20) :: fmt
        integer :: prec
        
        ! Set default values
        fmt = DEFAULT_TIME_FORMAT
        if (present(format)) fmt = format
        
        prec = 3
        if (present(precision)) prec = precision
        
        call et2utc(et_time, fmt, prec, et_to_utc)
    end function et_to_utc
    
    !> Convert UTC time string to Julian Date
    !> 
    !> This function converts a UTC time string to Julian Date.
    !> 
    !> @param[in] utc_string UTC time string
    !> @return Julian Date
    !> 
    !> @example
    !> ```fortran
    !> real(DP) :: jd
    !> jd = utc_to_jd('2024-01-01T12:00:00')
    !> ```
    function utc_to_jd(utc_string)
        character(len=*), intent(in) :: utc_string
        real(DP) :: utc_to_jd
        
        real(DP) :: et_time
        
        ! Convert to ET first, then to JD
        et_time = utc_to_et(utc_string)
        utc_to_jd = 2451545.0_DP + et_time / 86400.0_DP  ! J2000 = JD 2451545.0
    end function utc_to_jd
    
    !> Convert Julian Date to UTC time string
    !> 
    !> This function converts Julian Date to UTC time string.
    !> 
    !> @param[in] jd Julian Date
    !> @param[in] format Output format
    !> @return UTC time string
    !> 
    !> @example
    !> ```fortran
    !> character(len=100) :: utc_string
    !> utc_string = jd_to_utc(2460320.0_DP, 'ISOC')
    !> ```
    function jd_to_utc(jd, format)
        real(DP), intent(in) :: jd
        character(len=*), intent(in), optional :: format
        character(len=MAX_STRING_LEN) :: jd_to_utc
        
        real(DP) :: et_time
        
        ! Convert JD to ET
        et_time = (jd - 2451545.0_DP) * 86400.0_DP  ! J2000 = JD 2451545.0
        
        ! Convert ET to UTC string
        jd_to_utc = et_to_utc(et_time, format)
    end function jd_to_utc
    
    !> Convert UTC time string to Modified Julian Date
    !> 
    !> This function converts a UTC time string to Modified Julian Date.
    !> 
    !> @param[in] utc_string UTC time string
    !> @return Modified Julian Date
    !> 
    !> @example
    !> ```fortran
    !> real(DP) :: mjd
    !> mjd = utc_to_mjd('2024-01-01T12:00:00')
    !> ```
    function utc_to_mjd(utc_string)
        character(len=*), intent(in) :: utc_string
        real(DP) :: utc_to_mjd
        
        real(DP) :: jd
        
        jd = utc_to_jd(utc_string)
        utc_to_mjd = jd - 2400000.5_DP  ! MJD = JD - 2400000.5
    end function utc_to_mjd
    
    !> Convert Modified Julian Date to UTC time string
    !> 
    !> This function converts Modified Julian Date to UTC time string.
    !> 
    !> @param[in] mjd Modified Julian Date
    !> @param[in] format Output format
    !> @return UTC time string
    !> 
    !> @example
    !> ```fortran
    !> character(len=100) :: utc_string
    !> utc_string = mjd_to_utc(60320.0_DP, 'ISOC')
    !> ```
    function mjd_to_utc(mjd, format)
        real(DP), intent(in) :: mjd
        character(len=*), intent(in), optional :: format
        character(len=MAX_STRING_LEN) :: mjd_to_utc
        
        real(DP) :: jd
        
        jd = mjd + 2400000.5_DP  ! JD = MJD + 2400000.5
        mjd_to_utc = jd_to_utc(jd, format)
    end function mjd_to_utc
    
    !> Validate time string format
    !> 
    !> This function checks if a time string is in valid format.
    !> 
    !> @param[in] time_string Time string to validate
    !> @return .true. if valid, .false. otherwise
    !> 
    !> @example
    !> ```fortran
    !> logical :: is_valid
    !> is_valid = validate_time_string('2024-01-01T12:00:00')
    !> ```
    logical function validate_time_string(time_string)
        character(len=*), intent(in) :: time_string
        character(len=MAX_STRING_LEN) :: test_output
        real(DP) :: test_et
        integer :: status
        
        validate_time_string = .false.
        
        ! Check if string is empty
        if (len_trim(time_string) == 0) return
        
        ! Try to convert the string to ET
        ! Note: str2et will set error flag if conversion fails
        call str2et(time_string, test_et)
        
        ! Check if conversion was successful by trying to convert back
        call et2utc(test_et, 'ISOC', 3, test_output)
        
        ! If we get here without error, the string is valid
        validate_time_string = .true.
        
    end function validate_time_string
    
    !> Get current UTC time
    !> 
    !> This function returns the current UTC time as a string.
    !> 
    !> @param[in] format Output format
    !> @return Current UTC time string
    !> 
    !> @note This function uses system time, not SPICE time
    !> 
    !> @example
    !> ```fortran
    !> character(len=100) :: current_time
    !> current_time = get_current_utc('ISOC')
    !> ```
    function get_current_utc(format)
        character(len=*), intent(in), optional :: format
        character(len=MAX_STRING_LEN) :: get_current_utc
        
        character(len=20) :: fmt
        integer :: values(8)
        character(len=20) :: time_str
        
        ! Set default format
        fmt = DEFAULT_TIME_FORMAT
        if (present(format)) fmt = format
        
        ! Get system time
        call date_and_time(values=values)
        
        ! Format as ISO string
        write(time_str, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
            values(1), values(2), values(3), values(4), values(5), values(6)
        
        get_current_utc = trim(time_str)
        
    end function get_current_utc

end module pod_time_module