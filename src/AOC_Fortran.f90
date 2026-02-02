module AOC_Fortran
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, AOC_Fortran!"
  end subroutine say_hello
end module AOC_Fortran
