!> # CAT Frame Simple Module
!>
!> 轻量级的帧（frame）基础定义模块，仅提供数据类型与基础枚举，
!> 不直接依赖 SPICE 转换例程，方便被上层代码（测试/测量/对象模型）
!> 以最小依赖方式引用。
!>
!> ## Responsibilities（职责）
!> - 提供 `cat_frame_state` 轻量结构（位置、速度、帧名、历元、有效标记）
!> - 提供帧类别枚举（惯性/自转/地固/台站/轨道）
!> - 约定统一使用 SPICE 帧名字符串（如 `J2000`、`ITRF2000` 等）
!>
!> ## 与 `cat_frame_module` 的差异
!> - 本模块：仅含类型与基础枚举，无 SPICE 依赖；供轻量引用
!> - `cat_frame_module`：实现帧间转换与几何运算，集成 SPICE（`pxform`/`sxform`等）
!>
!> ## 命名与ID策略
!> - 推荐仅使用帧名字符串；如需帧ID，请在运行时通过 SPICE `namfrm(name,id)` 获取
!> - 不在程序中新增/维护本地帧ID常量，避免与内核不一致
!>
!> ## Dependencies
!> - `cat_global`: 基础数据类型与常量
!>
!> ## Version
!> - **Created**: 2025-09-09
!> - **Updated**: 2025-09-12 - 明确模块职责分离与SPICE帧名策略
!>
!> @note 若需要进行帧间转换，请使用 `cat_frame_module` 中的接口。
!> @warning 使用前请确保所用帧已在 SPICE FK 中正确定义并加载。

module cat_frame_simple_module
    use cat_global, only: DP, MAX_STRING_LEN
    implicit none
    
    private
    
    !> Frame type enumeration
    enum, bind(c)
        enumerator :: FRAME_TYPE_INERTIAL = 1
        enumerator :: FRAME_TYPE_ROTATING = 2
        enumerator :: FRAME_TYPE_BODY_FIXED = 3
        enumerator :: FRAME_TYPE_TOPOCENTRIC = 4
        enumerator :: FRAME_TYPE_ORBITAL = 5
    end enum
    
    !> Frame ID enumeration
    enum, bind(c)
        enumerator :: FRAME_ID_BCRS = 1      ! Barycentric Celestial Reference System
        enumerator :: FRAME_ID_GCRS = 2      ! Geocentric Celestial Reference System
        enumerator :: FRAME_ID_ITRF = 3      ! International Terrestrial Reference Frame
        enumerator :: FRAME_ID_EME2000 = 4   ! Earth Mean Equator and Equinox of J2000.0
        enumerator :: FRAME_ID_TEME = 5      ! True Equator Mean Equinox
        enumerator :: FRAME_ID_TOD = 6       ! True of Date
        enumerator :: FRAME_ID_TOPOCENTRIC = 7 ! Topocentric frame
        enumerator :: FRAME_ID_RTN = 8       ! Radial-Transverse-Normal
    end enum
    
    public :: FRAME_TYPE_INERTIAL, FRAME_TYPE_ROTATING, FRAME_TYPE_BODY_FIXED, &
              FRAME_TYPE_TOPOCENTRIC, FRAME_TYPE_ORBITAL, &
              FRAME_ID_BCRS, FRAME_ID_GCRS, FRAME_ID_ITRF, FRAME_ID_EME2000, &
              FRAME_ID_TEME, FRAME_ID_TOD, FRAME_ID_TOPOCENTRIC, FRAME_ID_RTN
    
    !> Frame state structure
    type, public :: cat_frame_state
        real(DP), dimension(3) :: position = 0.0_DP
        real(DP), dimension(3) :: velocity = 0.0_DP
        character(len=MAX_STRING_LEN) :: frame_name = ""
        real(DP) :: epoch = 0.0_DP
        logical :: valid = .false.
    contains
        procedure :: get_position
        procedure :: set_position
        procedure :: get_velocity
        procedure :: set_velocity
        procedure :: get_frame_name
        procedure :: set_frame_name
        procedure :: get_epoch
        procedure :: set_epoch
        procedure :: is_valid
        procedure :: set_valid
        procedure :: copy_state
        procedure :: clear_state
    end type cat_frame_state
    
contains
    
    ! ============================================================================
    ! Frame state procedures
    ! ============================================================================
    
    function get_position(self) result(position)
        class(cat_frame_state), intent(in) :: self
        real(DP), dimension(3) :: position
        position = self%position
    end function get_position
    
    subroutine set_position(self, position)
        class(cat_frame_state), intent(inout) :: self
        real(DP), dimension(3), intent(in) :: position
        self%position = position
    end subroutine set_position
    
    function get_velocity(self) result(velocity)
        class(cat_frame_state), intent(in) :: self
        real(DP), dimension(3) :: velocity
        velocity = self%velocity
    end function get_velocity
    
    subroutine set_velocity(self, velocity)
        class(cat_frame_state), intent(inout) :: self
        real(DP), dimension(3), intent(in) :: velocity
        self%velocity = velocity
    end subroutine set_velocity
    
    function get_frame_name(self) result(frame_name)
        class(cat_frame_state), intent(in) :: self
        character(len=MAX_STRING_LEN) :: frame_name
        frame_name = self%frame_name
    end function get_frame_name
    
    subroutine set_frame_name(self, frame_name)
        class(cat_frame_state), intent(inout) :: self
        character(len=*), intent(in) :: frame_name
        self%frame_name = frame_name
    end subroutine set_frame_name
    
    function get_epoch(self) result(epoch)
        class(cat_frame_state), intent(in) :: self
        real(DP) :: epoch
        epoch = self%epoch
    end function get_epoch
    
    subroutine set_epoch(self, epoch)
        class(cat_frame_state), intent(inout) :: self
        real(DP), intent(in) :: epoch
        self%epoch = epoch
    end subroutine set_epoch
    
    logical function is_valid(self)
        class(cat_frame_state), intent(in) :: self
        is_valid = self%valid
    end function is_valid
    
    subroutine set_valid(self, valid)
        class(cat_frame_state), intent(inout) :: self
        logical, intent(in) :: valid
        self%valid = valid
    end subroutine set_valid
    
    subroutine copy_state(self, other_state)
        class(cat_frame_state), intent(inout) :: self
        class(cat_frame_state), intent(in) :: other_state
        self%position = other_state%position
        self%velocity = other_state%velocity
        self%frame_name = other_state%frame_name
        self%epoch = other_state%epoch
        self%valid = other_state%valid
    end subroutine copy_state
    
    subroutine clear_state(self)
        class(cat_frame_state), intent(inout) :: self
        self%position = 0.0_DP
        self%velocity = 0.0_DP
        self%frame_name = ""
        self%epoch = 0.0_DP
        self%valid = .false.
    end subroutine clear_state

end module cat_frame_simple_module
