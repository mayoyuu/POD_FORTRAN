!> # CAT Object Base Module
!> 
!> This module provides the abstract base class and interfaces for all CAT objects
!> in the object-oriented architecture.
!> 
!> ## Features
!> 
!> - **Abstract Base Class**: Common interface for all CAT objects
!> - **Type Extensions**: Support for inheritance and polymorphism
!> - **Generic Interfaces**: Type-safe procedure interfaces
!> - **Object Management**: Lifecycle management for objects
!> 
!> ## Dependencies
!> 
!> - `pod_global`: For data types and constants
!> 
!> ## Author
!> 
!> Zhao Yuhui (PMO, zhaoyuhui@pmo.ac.cn)
!> 
!> ## Version
!> 
!> - **Created**: 2025-09-09
!> - **Ref**: [待添加参考文献]
!> 
!> @note This module provides the foundation for the object-oriented CAT system
!> @warning All CAT objects must inherit from pod_object

module pod_object_base
    use pod_global, only: DP, MAX_STRING_LEN
    implicit none
    
    private
    
    !> Abstract base class for all CAT objects
    type, abstract, public :: pod_object
        character(len=MAX_STRING_LEN) :: name = ""
        character(len=MAX_STRING_LEN) :: id = ""
        logical :: initialized = .false.
        integer :: error_code = 0
        character(len=MAX_STRING_LEN) :: error_message = ""
    contains
        ! Abstract procedures that must be implemented by derived types
        procedure(init_interface), deferred :: initialize
        procedure(cleanup_interface), deferred :: cleanup
        procedure(is_initialized_interface), deferred :: is_initialized
        procedure(get_name_interface), deferred :: get_name
        procedure(set_name_interface), deferred :: set_name
        
        ! Common procedures with default implementations
        procedure :: get_id
        procedure :: set_id
        procedure :: get_error_code
        procedure :: set_error
        procedure :: clear_error
        procedure :: has_error
        procedure :: get_error_message
    end type pod_object
    
    ! Abstract interfaces for deferred procedures
    abstract interface
        !> Initialize the object
        subroutine init_interface(self, config_file)
            import :: pod_object, MAX_STRING_LEN
            class(pod_object), intent(inout) :: self
            character(len=*), intent(in), optional :: config_file
        end subroutine init_interface
        
        !> Cleanup the object
        subroutine cleanup_interface(self)
            import :: pod_object
            class(pod_object), intent(inout) :: self
        end subroutine cleanup_interface
        
        !> Check if object is initialized
        logical function is_initialized_interface(self)
            import :: pod_object
            class(pod_object), intent(in) :: self
        end function is_initialized_interface
        
        !> Get object name
        function get_name_interface(self) result(name)
            import :: pod_object, MAX_STRING_LEN
            class(pod_object), intent(in) :: self
            character(len=MAX_STRING_LEN) :: name
        end function get_name_interface
        
        !> Set object name
        subroutine set_name_interface(self, name)
            import :: pod_object, MAX_STRING_LEN
            class(pod_object), intent(inout) :: self
            character(len=*), intent(in) :: name
        end subroutine set_name_interface
    end interface
    
    !> Object pointer type for polymorphic objects
    type, public :: pod_object_ptr
        class(pod_object), pointer :: obj => null()
    contains
        procedure :: is_associated
        procedure :: get_object
        procedure :: set_object
        procedure :: destroy
    end type pod_object_ptr
    
    !> Object factory interface
    type, abstract, public :: pod_factory
        character(len=MAX_STRING_LEN) :: factory_name = ""
    contains
        procedure(create_object_interface), deferred :: create_object
        procedure :: get_factory_name
        procedure :: set_factory_name
    end type pod_factory
    
    abstract interface
        !> Create an object of specified type
        function create_object_interface(self, object_type, name) result(obj)
            import :: pod_factory, pod_object_ptr, MAX_STRING_LEN
            class(pod_factory), intent(in) :: self
            character(len=*), intent(in) :: object_type
            character(len=*), intent(in), optional :: name
            type(pod_object_ptr) :: obj
        end function create_object_interface
    end interface
    
contains
    
    !> Get object ID
    function get_id(self) result(id)
        class(pod_object), intent(in) :: self
        character(len=MAX_STRING_LEN) :: id
        id = self%id
    end function get_id
    
    !> Set object ID
    subroutine set_id(self, id)
        class(pod_object), intent(inout) :: self
        character(len=*), intent(in) :: id
        self%id = id
    end subroutine set_id
    
    !> Get error code
    function get_error_code(self) result(error_code)
        class(pod_object), intent(in) :: self
        integer :: error_code
        error_code = self%error_code
    end function get_error_code
    
    !> Set error
    subroutine set_error(self, error_code, error_message)
        class(pod_object), intent(inout) :: self
        integer, intent(in) :: error_code
        character(len=*), intent(in) :: error_message
        self%error_code = error_code
        self%error_message = error_message
    end subroutine set_error
    
    !> Clear error
    subroutine clear_error(self)
        class(pod_object), intent(inout) :: self
        self%error_code = 0
        self%error_message = ""
    end subroutine clear_error
    
    !> Check if object has error
    logical function has_error(self)
        class(pod_object), intent(in) :: self
        has_error = (self%error_code /= 0)
    end function has_error
    
    !> Get error message
    function get_error_message(self) result(error_message)
        class(pod_object), intent(in) :: self
        character(len=MAX_STRING_LEN) :: error_message
        error_message = self%error_message
    end function get_error_message
    
    !> Check if object pointer is associated
    logical function is_associated(self)
        class(pod_object_ptr), intent(in) :: self
        is_associated = associated(self%obj)
    end function is_associated
    
    !> Get object from pointer
    function get_object(self) result(obj)
        class(pod_object_ptr), intent(in) :: self
        class(pod_object), pointer :: obj
        obj => self%obj
    end function get_object
    
    !> Set object in pointer
    subroutine set_object(self, obj)
        class(pod_object_ptr), intent(inout) :: self
        class(pod_object), intent(in), target :: obj
        self%obj => obj
    end subroutine set_object
    
    !> Destroy object
    subroutine destroy(self)
        class(pod_object_ptr), intent(inout) :: self
        if (associated(self%obj)) then
            call self%obj%cleanup()
            deallocate(self%obj)
            self%obj => null()
        end if
    end subroutine destroy
    
    !> Get factory name
    function get_factory_name(self) result(name)
        class(pod_factory), intent(in) :: self
        character(len=MAX_STRING_LEN) :: name
        name = self%factory_name
    end function get_factory_name
    
    !> Set factory name
    subroutine set_factory_name(self, name)
        class(pod_factory), intent(inout) :: self
        character(len=*), intent(in) :: name
        self%factory_name = name
    end subroutine set_factory_name

end module pod_object_base
