module vector
  implicit none
  type,public :: vector3
     sequence
     real(kind=8) :: vec(3)
  end type vector3

contains

  subroutine vector3_read(FILE, v3)
    integer, intent(IN) :: FILE
    type(vector3), intent(OUT) :: v3
    read(FILE,*) v3%vec(1), v3%vec(2), v3%vec(3)
  end subroutine vector3_read

  function vector3_inner_product(v1,v2)
    real(kind=8)             :: vector3_inner_product
    type(vector3),intent(in) :: v1,v2
    !local variables
    real(kind=8)             :: sum
    integer                  :: k
    sum = 0
    do k=1,3
       sum = sum + v1%vec(k) * v2%vec(k)
    enddo
    vector3_inner_product = sum
  end function vector3_inner_product

end module vector
